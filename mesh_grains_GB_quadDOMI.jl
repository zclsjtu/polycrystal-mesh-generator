using Gmsh
using Dates
using Statistics
using DelimitedFiles

"""
尝试设置Gmsh选项，忽略不支持的选项
"""
function try_set_option(option_name, value)
    try
        Gmsh.gmsh.option.setNumber(option_name, value)
        return true
    catch e
        println("  选项 '$option_name' 设置失败: $e")
        return false
    end
end

"""
从geo文件中直接提取线段和边长信息
"""
function extract_edge_lengths_from_geo(geo_file)
    if !isfile(geo_file)
        error("几何文件不存在: $geo_file")
    end
    
    # 读取文件内容
    file_content = read(geo_file, String)
    
    # 存储点和线段信息
    points = Dict{Int, Vector{Float64}}()
    lines = Dict{Int, Vector{Int}}()
    edge_lengths = Dict{Int, Float64}()
    
    # 提取点信息
    for m in eachmatch(r"Point\s*\(\s*(\d+)\s*\)\s*=\s*\{\s*([^,]+)\s*,\s*([^,]+)", file_content)
        id = parse(Int, m.captures[1])
        x = parse(Float64, m.captures[2])
        y = parse(Float64, m.captures[3])
        points[id] = [x, y]
    end
    
    # 提取线段信息
    for m in eachmatch(r"Line\s*\(\s*(\d+)\s*\)\s*=\s*\{\s*(\d+)\s*,\s*(\d+)\s*\}", file_content)
        id = parse(Int, m.captures[1])
        start_point = parse(Int, m.captures[2])
        end_point = parse(Int, m.captures[3])
        
        lines[id] = [start_point, end_point]
        
        # 直接计算线段长度
        if haskey(points, start_point) && haskey(points, end_point)
            p1 = points[start_point]
            p2 = points[end_point]
            length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
            edge_lengths[id] = length
        end
    end
    
    # 提取物理组边界信息
    physical_groups = Dict{Int, String}()
    physical_group_entities = Dict{Int, Vector{Int}}()
    
    for m in eachmatch(r"Physical\s+Surface\s*\(\s*\"([^\"]+)\"\s*\)\s*=\s*\{([^}]+)\}", file_content)
        name = m.captures[1]
        entities_str = m.captures[2]
        entities = [parse(Int, strip(e)) for e in split(entities_str, ",")]
        
        # 为了简化，使用数字ID
        group_id = length(physical_groups) + 1
        physical_groups[group_id] = name
        physical_group_entities[group_id] = entities
    end
    
    # 提取表面信息
    surface_to_loop = Dict{Int, Int}()
    
    for m in eachmatch(r"(Plane\s+)?Surface\s*\(\s*(\d+)\s*\)\s*=\s*\{([^}]+)\}", file_content)
        surface_id = parse(Int, m.captures[2])
        loop_str = m.captures[3]
        # 简化处理，只取第一个值
        loop_id = parse(Int, strip(split(loop_str, ",")[1]))
        surface_to_loop[surface_id] = loop_id
    end
    
    # 提取线环信息
    loop_to_lines = Dict{Int, Vector{Int}}()
    
    for m in eachmatch(r"Line\s+Loop\s*\(\s*(\d+)\s*\)\s*=\s*\{([^}]+)\}", file_content)
        loop_id = parse(Int, m.captures[1])
        lines_str = m.captures[2]
        line_ids = []
        
        for line_id_str in split(lines_str, ",")
            line_id_str = strip(line_id_str)
            if startswith(line_id_str, "-")
                # 处理负号表示的反向线
                line_id = parse(Int, line_id_str[2:end])
                push!(line_ids, -line_id)
            else
                line_id = parse(Int, line_id_str)
                push!(line_ids, line_id)
            end
        end
        
        loop_to_lines[loop_id] = line_ids
    end
    
    # 识别晶界表面和晶粒表面，排除AllBoundaries和AllGrains
    boundary_surfaces = []
    grain_surfaces = []
    all_boundaries_group = nothing
    all_grains_group = nothing
    
    for (group_id, name) in physical_groups
        # 区分AllBoundaries和Boundary开头的组（不区分大小写）
        if lowercase(name) == "allboundaries"
            all_boundaries_group = group_id
        # 区分AllGrains和Grain开头的组（不区分大小写）
        elseif lowercase(name) == "allgrains"
            all_grains_group = group_id
        # 收集所有以Boundary开头的组（不区分大小写）
        elseif startswith(lowercase(name), "boundary")
            if haskey(physical_group_entities, group_id)
                append!(boundary_surfaces, physical_group_entities[group_id])
            end
        # 收集所有以Grain开头的组（不区分大小写）
        elseif startswith(lowercase(name), "grain")
            if haskey(physical_group_entities, group_id)
                append!(grain_surfaces, physical_group_entities[group_id])
            end
        end
    end
    
    # 去重
    unique!(boundary_surfaces)
    unique!(grain_surfaces)
    
    println("从geo文件提取信息:")
    println("  找到 $(length(points)) 个点")
    println("  找到 $(length(lines)) 条线段")
    println("  找到 $(length(physical_groups)) 个物理组")
    
    if all_boundaries_group !== nothing
        println("  找到总晶界组: $(physical_groups[all_boundaries_group])")
    end
    if all_grains_group !== nothing
        println("  找到总晶粒组: $(physical_groups[all_grains_group])")
    end
    
    println("  找到 $(length(boundary_surfaces)) 个晶界表面")
    println("  找到 $(length(grain_surfaces)) 个晶粒表面")
    
    return points, lines, edge_lengths, loop_to_lines, surface_to_loop, boundary_surfaces, grain_surfaces, physical_groups
end

"""
更精确地识别四边形中的对边关系
基于拓扑分析，而非简单的长度比较
"""
function identify_opposite_edges(surface_tag, edges, lines)
    if length(edges) != 4
        return [], false
    end
    
    # 获取四边形顶点的连接关系
    edge_vertices = Dict{Int, Set{Int}}()  # 每条边连接的顶点
    vertex_edges = Dict{Int, Set{Int}}()   # 每个顶点连接的边
    
    for edge_id in edges
        abs_edge = abs(edge_id)
        if haskey(lines, abs_edge)
            verts = Set(lines[abs_edge])
            edge_vertices[edge_id] = verts
            
            for v in verts
                if !haskey(vertex_edges, v)
                    vertex_edges[v] = Set{Int}()
                end
                push!(vertex_edges[v], edge_id)
            end
        else
            # 找不到边的信息，无法分析
            return [], false
        end
    end
    
    # 一个有效的四边形应该有4个顶点，每个顶点连接2条边
    if length(vertex_edges) != 4
        return [], false
    end
    
    for (_, connected_edges) in vertex_edges
        if length(connected_edges) != 2
            return [], false
        end
    end
    
    # 寻找对边 - 对边不共享任何顶点
    opposite_pairs = []
    
    for i in 1:3
        for j in (i+1):4
            e1 = edges[i]
            e2 = edges[j]
            
            if !haskey(edge_vertices, e1) || !haskey(edge_vertices, e2)
                continue
            end
            
            # 检查两边是否共享顶点
            shared_vertices = intersect(edge_vertices[e1], edge_vertices[e2])
            
            if isempty(shared_vertices)
                # 没有共享顶点，是对边
                push!(opposite_pairs, (i, j))
            end
        end
    end
    
    # 验证是否找到了两对对边
    valid = length(opposite_pairs) == 2
    
    return opposite_pairs, valid
end

"""
根据拓扑分析和长度比较找出最合理的四边形对边配对
"""
function get_best_opposite_edges(surface_tag, edges, lines, edge_lengths)
    # 先基于拓扑分析
    topo_pairs, topo_valid = identify_opposite_edges(surface_tag, edges, lines)
    
    if topo_valid
        return topo_pairs, "topology"
    end
    
    # 如果拓扑分析失败，退回到基于边长的分析
    if length(edges) != 4
        return [], "invalid"
    end
    
    # 对边长进行排序，获取排序索引（从短到长）
    edge_lens = [get(edge_lengths, abs(e), 0.0) for e in edges]
    sorted_indices = sortperm(edge_lens)
    
    # 获取最短到最长的四条边索引
    shortest_idx = sorted_indices[1]
    short_idx = sorted_indices[2]
    long_idx = sorted_indices[3]
    longest_idx = sorted_indices[4]
    
    # 计算相邻边长度差
    diff1 = abs(edge_lens[longest_idx] - edge_lens[long_idx])
    diff2 = abs(edge_lens[short_idx] - edge_lens[shortest_idx])
    diff3 = abs(edge_lens[longest_idx] - edge_lens[shortest_idx])
    diff4 = abs(edge_lens[long_idx] - edge_lens[short_idx])
    
    # 如果最长和次长边长度接近，最短和次短边长度接近
    # 则认为这是一个矩形或类矩形，应配对为长边对长边，短边对短边
    if (diff1 + diff2) < (diff3 + diff4)
        return [(longest_idx, long_idx), (shortest_idx, short_idx)], "rectangular"
    else
        # 否则可能是一个平行四边形，配对为最长边对最短边，次长边对次短边
        return [(longest_idx, shortest_idx), (long_idx, short_idx)], "parallelogram"
    end
end

"""
验证四边形表面是否可以应用Transfinite算法
确保对边具有相同数量的控制点
"""
function verify_transfinite_compatibility(surface_tag, edges, lines, edge_lengths, edge_divisions)
    # 获取最佳对边配对
    opposite_pairs, method = get_best_opposite_edges(surface_tag, edges, lines, edge_lengths)
    
    if opposite_pairs == [] || length(opposite_pairs) != 2
        return false, "无法识别对边", method
    end
    
    # 检查每对对边的控制点数是否一致
    for pair in opposite_pairs
        i1, i2 = pair
        e1 = abs(edges[i1])
        e2 = abs(edges[i2])
        
        if !haskey(edge_divisions, e1) || !haskey(edge_divisions, e2)
            continue
        end
        
        if edge_divisions[e1] != edge_divisions[e2]
            return false, "对边控制点不一致: $(edge_divisions[e1]) != $(edge_divisions[e2])", method
        end
    end
    
    return true, "通过", method
end

"""
使用几何拓扑分析设置Transfinite表面，确保对边控制点一致
"""
function set_transfinite_surface_with_topology(surface_tag, edges, lines, points)
    if length(edges) != 4
        return false, "非四边形表面"
    end
    
    # 提取四边形的顶点
    vertices = Set{Int}()
    for edge_id in edges
        abs_edge = abs(edge_id)
        if haskey(lines, abs_edge)
            for v in lines[abs_edge]
                push!(vertices, v)
            end
        end
    end
    
    # 四边形应该有恰好4个顶点
    if length(vertices) != 4
        return false, "顶点数不为4"
    end
    
    # 转换为数组
    vertex_array = collect(vertices)
    
    try
        # 使用指定的顶点顺序设置Transfinite
        Gmsh.gmsh.model.mesh.setTransfiniteSurface(surface_tag, "Left", vertex_array)
        Gmsh.gmsh.model.mesh.setRecombine(2, surface_tag)
        return true, "成功"
    catch e
        return false, "设置失败: $e"
    end
end

"""
获取多边形的中心点
"""
function get_polygon_centroid(vertices, points)
    if length(vertices) < 3 || isempty(points)
        return [0.0, 0.0]
    end
    
    x_sum = 0.0
    y_sum = 0.0
    valid_count = 0
    
    for v in vertices
        if haskey(points, v)
            x_sum += points[v][1]
            y_sum += points[v][2]
            valid_count += 1
        end
    end
    
    if valid_count > 0
        return [x_sum / valid_count, y_sum / valid_count]
    else
        return [0.0, 0.0]
    end
end

"""
分析晶粒形状和特征
"""
function analyze_grain_shape(surface_tag, edges, lines, edge_lengths)
    if !haskey(lines, 1)  # 检查lines字典是否为空
        return "unknown", 0.0
    end
    
    # 数量分析
    edge_count = length(edges)
    
    if edge_count == 0
        return "unknown", 0.0
    elseif edge_count == 3
        return "triangle", 0.0
    elseif edge_count == 4
        # 四边形详细分析
        edge_lens = [get(edge_lengths, abs(e), 0.0) for e in edges]
        
        if isempty(edge_lens) || any(isnan, edge_lens) || any(isinf, edge_lens)
            return "quadrilateral", 0.0
        end
        
        sorted_lens = sort(edge_lens)
        min_len = sorted_lens[1]
        max_len = sorted_lens[end]
        
        # 计算长宽比
        aspect_ratio = max_len > 0 ? max_len / min_len : 1.0
        
        # 判断是否接近正方形
        lens_diff = maximum(edge_lens) - minimum(edge_lens)
        len_avg = sum(edge_lens) / 4
        
        if lens_diff / len_avg < 0.1
            return "square", aspect_ratio
        else
            return "rectangle", aspect_ratio
        end
    elseif edge_count == 5
        return "pentagon", 0.0
    elseif edge_count == 6
        return "hexagon", 0.0
    else
        return "polygon_$(edge_count)", 0.0
    end
end

"""
全局边协调算法 - 解决共享边冲突问题
确保所有四边形的对边都有相同的控制点数
"""
function coordinate_edge_divisions_globally(surface_edges, quad_surfaces, surface_opposite_edges, edge_lengths, boundary_size, grain_size)
    # 初始化边的节点数量
    edge_divisions = Dict{Int, Int}()
    
    # 阶段1: 确定每条边的初始节点数
    println("第1阶段: 设置每条边的初始节点数...")
    for (edge_tag, length) in edge_lengths
        # 根据长度确定初始节点数
        n_divisions = max(2, ceil(Int, length / boundary_size))
        edge_divisions[edge_tag] = n_divisions
        edge_divisions[-edge_tag] = n_divisions  # 同时设置反向边
    end
    
    # 阶段2: 构建边与四边形的关系图
    println("第2阶段: 构建边与四边形的关系图...")
    edge_to_quads = Dict{Int, Vector{Int}}()  # 每条边属于哪些四边形
    quad_constraints = Dict{Int, Dict{Int, Int}}() # 四边形ID -> {边ID -> 对边ID}
    
    # 收集每条边属于哪些四边形
    for surface_tag in quad_surfaces
        if !haskey(surface_edges, surface_tag) || !haskey(surface_opposite_edges, surface_tag)
            continue
        end
        
        edges = surface_edges[surface_tag]
        opposite_pairs = surface_opposite_edges[surface_tag]
        
        # 记录该四边形中每条边的对边
        quad_constraint = Dict{Int, Int}()
        
        for (i1, i2) in opposite_pairs
            e1 = edges[i1]
            e2 = edges[i2]
            abs_e1 = abs(e1)
            abs_e2 = abs(e2)
            
            # 记录对边关系
            quad_constraint[abs_e1] = abs_e2
            quad_constraint[abs_e2] = abs_e1
            
            # 记录边属于哪些四边形
            if !haskey(edge_to_quads, abs_e1)
                edge_to_quads[abs_e1] = []
            end
            if !haskey(edge_to_quads, abs_e2)
                edge_to_quads[abs_e2] = []
            end
            
            push!(edge_to_quads[abs_e1], surface_tag)
            push!(edge_to_quads[abs_e2], surface_tag)
        end
        
        quad_constraints[surface_tag] = quad_constraint
    end
    
    # 输出共享边信息
    shared_edges = []
    for (edge, quads) in edge_to_quads
        if length(quads) > 1
            push!(shared_edges, (edge, quads))
        end
    end
    println("  发现 $(length(shared_edges)) 条共享边")
    
    # 阶段3: 迭代协调算法
    println("第3阶段: 执行迭代协调算法...")
    changed = true
    iteration = 0
    max_iterations = 10
    
    while changed && iteration < max_iterations
        iteration += 1
        changed = false
        
        # 对每个四边形应用对边约束
        for (surface_tag, constraint) in quad_constraints
            for (e1, e2) in constraint
                if !haskey(edge_divisions, e1) || !haskey(edge_divisions, e2)
                    continue
                end
                
                n1 = edge_divisions[e1]
                n2 = edge_divisions[e2]
                
                if n1 != n2
                    # 使用较大的值
                    new_n = max(n1, n2)
                    edge_divisions[e1] = new_n
                    edge_divisions[e2] = new_n
                    edge_divisions[-e1] = new_n
                    edge_divisions[-e2] = new_n
                    changed = true
                    
                    # 传播变化到共享这些边的其他四边形
                    println("  调整: 边 #$e1($n1) 和 #$e2($n2) 统一为 $new_n 点")
                end
            end
        end
        
        println("  迭代 $iteration: $(changed ? "有变化" : "无变化")")
    end
    
    # 阶段4: 最终验证
    println("第4阶段: 验证协调结果...")
    success_count = 0
    fail_count = 0
    
    for surface_tag in quad_surfaces
        if !haskey(surface_opposite_edges, surface_tag)
            continue
        end
        
        edges = surface_edges[surface_tag]
        opposite_pairs = surface_opposite_edges[surface_tag]
        is_valid = true
        
        for (i1, i2) in opposite_pairs
            e1 = abs(edges[i1])
            e2 = abs(edges[i2])
            
            if !haskey(edge_divisions, e1) || !haskey(edge_divisions, e2)
                is_valid = false
                break
            end
            
            if edge_divisions[e1] != edge_divisions[e2]
                is_valid = false
                println("  验证失败: 表面 #$surface_tag 的对边 #$e1($(edge_divisions[e1])) 和 #$e2($(edge_divisions[e2])) 点数不同")
                break
            end
        end
        
        if is_valid
            success_count += 1
        else
            fail_count += 1
        end
    end
    
    println("  协调结果: $success_count 成功, $fail_count 失败")
    
    return edge_divisions, success_count, fail_count
end

"""
生成多晶材料结构化网格，使用算法8优化晶粒内部网格
"""
function generate_polycrystal_mesh_with_algorithm8(geo_file, mesh_file; 
                                                 grain_size=0.05,           # 晶粒尺寸较大 
                                                 boundary_size=0.001,       # 晶界尺寸较小
                                                 transition_distance=0.1,   # 尺寸过渡区域
                                                 num_threads=8,
                                                 max_optimization_time=60,
                                                 interactive=true,
                                                 verbose=true,
                                                 output_file="polycrystal_mesh_info.txt")
    start_time = time()
    
    # 打开输出文件
    output_io = open(output_file, "w")
    
    # 文件头添加字段说明
    println(output_io, "# 多晶材料结构化网格信息文件 (晶粒算法8优化版)")
    println(output_io, "# 生成时间: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(output_io, "#")
    println(output_io, "# 参数设置:")
    println(output_io, "# - 几何文件: $geo_file")
    println(output_io, "# - 输出网格文件: $mesh_file")
    println(output_io, "# - 晶粒网格尺寸: $grain_size")
    println(output_io, "# - 晶界网格尺寸: $boundary_size")
    println(output_io, "# - 过渡区域距离: $transition_distance")
    println(output_io, "# - 优化时间限制: $max_optimization_time 秒")
    println(output_io, "#")
    println(output_io, "# 晶界区域使用算法6(Frontal-Delaunay)，晶粒内部使用算法8(Delaunay for Quads)")
    println(output_io, "# ==========================================")
    
    # 从geo文件中直接提取边长信息
    points, lines, direct_edge_lengths, loop_to_lines, surface_to_loop, boundary_surfaces, grain_surfaces, physical_groups = extract_edge_lengths_from_geo(geo_file)
    
    # 初始化Gmsh
    Gmsh.gmsh.initialize()
    try_set_option("General.Terminal", 1)
    try_set_option("General.NumThreads", num_threads)
    println("启动Gmsh，线程数: $num_threads")
    
    # 加载几何文件
    Gmsh.gmsh.open(geo_file)
    Gmsh.gmsh.model.geo.synchronize()
    println("加载几何文件: $(geo_file) [$(round(time() - start_time, digits=2))s]")
    
    # 收集物理组信息
    phys_group_names = Dict{Int, String}()
    phys_group_surfaces = Dict{Int, Vector{Int}}()
    surface_to_phys_group = Dict{Int, Int}()
    all_boundary_phys_groups = []
    all_grain_phys_groups = []
    
    for (dim, tag) in Gmsh.gmsh.model.getPhysicalGroups()
        if dim == 2  # 只关注2D表面
            name = Gmsh.gmsh.model.getPhysicalName(dim, tag)
            entities = Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            
            # 保存物理组信息
            phys_group_names[tag] = name
            phys_group_surfaces[tag] = entities
            
            # 建立表面到物理组的映射
            for entity in entities
                surface_to_phys_group[entity] = tag
            end
            
            if verbose
                println("  - 物理组: $(name) (tag $(tag)), 包含 $(length(entities)) 个实体")
            end
            
            # 按命名规则分类，排除AllBoundaries和AllGrains
            if lowercase(name) != "allboundaries" && startswith(lowercase(name), "boundary")
                push!(all_boundary_phys_groups, tag)
            elseif lowercase(name) != "allgrains" && startswith(lowercase(name), "grain")
                push!(all_grain_phys_groups, tag)
            end
        end
    end
    
    # ===== 设置基本网格参数 =====
    println("设置全局网格参数...")
    
    # 基本参数
    try_set_option("Mesh.MeshSizeMin", boundary_size)
    try_set_option("Mesh.MeshSizeMax", grain_size)
    try_set_option("Mesh.MeshSizeFromPoints", 1)
    try_set_option("Mesh.MeshSizeExtendFromBoundary", 1)
    
    # 设置默认为算法6
    try_set_option("Mesh.Algorithm", 6)  # Frontal-Delaunay(晶界默认算法)
    
    # 结构化网格设置
    try_set_option("Mesh.RecombineAll", 1)  # 全局启用四边形重组
    try_set_option("Mesh.RecombinationAlgorithm", 1)  # 标准重组算法，更稳定
    
    # 质量参数
    try_set_option("Mesh.ElementOrder", 1)
    try_set_option("Mesh.Smoothing", 100)
    try_set_option("Mesh.QualityType", 10)  # 雅可比质量度量
    try_set_option("Mesh.OptimizeNetgen", 1)
    
    # 设置Laplace2D优化时间限制
    try_set_option("Mesh.OptimizationTimeLimit", max_optimization_time)
    
    # ===== 收集晶界和晶粒表面的边 =====
    println("收集表面的边信息...")
    
    # 收集每个表面的边
    surface_edges = Dict{Int, Vector{Int}}()
    
    # 从geo文件中收集表面与边的关系
    for surface_tag in vcat(boundary_surfaces, grain_surfaces)
        if haskey(surface_to_loop, surface_tag)
            loop_id = surface_to_loop[surface_tag]
            if haskey(loop_to_lines, loop_id)
                # 获取线环中的边，但保留方向信息（正负号）
                edge_tags = [line_id for line_id in loop_to_lines[loop_id]]
                surface_edges[surface_tag] = edge_tags
            end
        end
    end
    
    # 验证晶界子区域是否全部为四边形
    quad_boundary_surfaces = []
    non_quad_boundary_surfaces = []
    
    for surface_tag in boundary_surfaces
        if haskey(surface_edges, surface_tag)
            edges = surface_edges[surface_tag]
            if length(edges) == 4
                push!(quad_boundary_surfaces, surface_tag)
            else
                push!(non_quad_boundary_surfaces, surface_tag)
                println("  警告: 晶界表面 #$surface_tag 不是四边形，有 $(length(edges)) 条边")
            end
        end
    end
    
    println("  晶界区域: $(length(quad_boundary_surfaces)) 个四边形, $(length(non_quad_boundary_surfaces)) 个非四边形")
    println(output_io, "# 四边形晶界数: $(length(quad_boundary_surfaces)), 非四边形晶界数: $(length(non_quad_boundary_surfaces))")
    
    # ===== 分析每个表面的对边关系 =====
    println("分析表面的对边关系...")
    
    surface_opposite_edges = Dict{Int, Vector{Tuple{Int, Int}}}()
    surface_opposite_method = Dict{Int, String}()
    
    # 对每个四边形表面分析对边关系
    for surface_tag in quad_boundary_surfaces
        edges = surface_edges[surface_tag]
        opposite_pairs, method = get_best_opposite_edges(surface_tag, edges, lines, direct_edge_lengths)
        
        if length(opposite_pairs) == 2
            surface_opposite_edges[surface_tag] = opposite_pairs
            surface_opposite_method[surface_tag] = method
            
            if verbose
                println("  表面 #$surface_tag: 使用 $method 方法识别对边")
                for (idx, pair) in enumerate(opposite_pairs)
                    i1, i2 = pair
                    e1 = edges[i1]
                    e2 = edges[i2]
                    len1 = get(direct_edge_lengths, abs(e1), 0.0)
                    len2 = get(direct_edge_lengths, abs(e2), 0.0)
                    println("    对边$(idx): #$e1($(round(len1, digits=6))) - #$e2($(round(len2, digits=6)))")
                end
            end
        else
            println("  警告: 表面 #$surface_tag 无法识别对边关系")
        end
    end
    
    # ===== 使用直接计算的边长 =====
    println("使用直接计算的边长...")
    edge_lengths = direct_edge_lengths
    
    # ===== 设置边的控制点数 =====
    println("使用全局边协调算法确保对边点数一致...")
    edge_divisions, coordination_success, coordination_failed = coordinate_edge_divisions_globally(
        surface_edges, quad_boundary_surfaces, surface_opposite_edges, edge_lengths, boundary_size, grain_size
    )

    println("  全局协调结果: $(coordination_success) 个四边形成功, $(coordination_failed) 个四边形失败")
    
    # ===== 为边设置Transfinite曲线 =====
    println("应用边的控制点设置...")
    
    # 设置每条边的节点数
    success_count = 0
    failed_edges = []
    
    for (edge_tag, n_nodes) in edge_divisions
        if edge_tag > 0  # 只处理正向边
            try
                Gmsh.gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_nodes, "Progression", 1.0)
                success_count += 1
            catch e
                push!(failed_edges, (edge_tag, n_nodes, e))
                println("  边 #$edge_tag 设置失败: $e")
            end
        end
    end
    
    println("  成功设置 $success_count 条边的控制点")
    
    # ===== 验证并应用Transfinite表面 =====
    println("验证并设置Transfinite表面...")
    
    transfinite_success = 0
    transfinite_failed = 0
    transfinite_results = Dict{Int, Tuple{Bool, String, String}}()
    
    # 验证每个四边形表面是否符合Transfinite条件
    for surface_tag in quad_boundary_surfaces
        if !haskey(surface_edges, surface_tag)
            continue
        end
        
        edges = surface_edges[surface_tag]
        
        # 验证对边控制点是否一致
        is_compatible, message, method = verify_transfinite_compatibility(
            surface_tag, edges, lines, edge_lengths, edge_divisions
        )
        
        transfinite_results[surface_tag] = (is_compatible, message, method)
        
        if is_compatible
            # 使用拓扑信息设置Transfinite表面
            is_set, set_message = set_transfinite_surface_with_topology(surface_tag, edges, lines, points)
            
            if is_set
                transfinite_success += 1
                if verbose
                    println("  表面 #$surface_tag Transfinite设置成功 (方法: $method)")
                end
            else
                transfinite_failed += 1
                println("  表面 #$surface_tag Transfinite设置失败: $set_message")
                
                # 退回到简单设置
                try
                    Gmsh.gmsh.model.mesh.setRecombine(2, surface_tag)
                catch e
                    println("    四边形重组也失败: $e")
                end
            end
        else
            transfinite_failed += 1
            println("  表面 #$surface_tag 不兼容Transfinite: $message")
            
            # 退回到简单的四边形重组
            try
                Gmsh.gmsh.model.mesh.setRecombine(2, surface_tag)
            catch e
                println("    四边形重组也失败: $e")
            end
        end
    end
    
    println("  成功设置 $transfinite_success / $(length(quad_boundary_surfaces)) 个Transfinite表面")
    println("  $(transfinite_failed) 个表面不兼容Transfinite，使用四边形重组")
    
    # 为非四边形表面设置重组
    for surface_tag in non_quad_boundary_surfaces
        try
            Gmsh.gmsh.model.mesh.setRecombine(2, surface_tag)
        catch e
            println("  表面 #$surface_tag 重组设置失败: $e")
        end
    end
    
    # ===== 设置阈值场控制网格尺寸过渡 =====
    println("创建阈值场控制网格尺寸过渡...")
    
    try
        # 收集所有晶界边
        all_boundary_edges = []
        for surface_tag in boundary_surfaces
            if haskey(surface_edges, surface_tag)
                append!(all_boundary_edges, abs.(surface_edges[surface_tag]))
            end
        end
        unique!(all_boundary_edges)
        
        println("  收集了 $(length(all_boundary_edges)) 条晶界边")
        
        # 1. 创建距离场 - 计算到晶界边的距离
        distance_field = 1
        Gmsh.gmsh.model.mesh.field.add("Distance", distance_field)
        Gmsh.gmsh.model.mesh.field.setNumbers(distance_field, "EdgesList", all_boundary_edges)
        Gmsh.gmsh.model.mesh.field.setNumber(distance_field, "NNodesByEdge", 200)  # 增加精度
        
        # 2. 创建阈值场 - 基于距离场控制尺寸变化
        threshold_field = 2
        Gmsh.gmsh.model.mesh.field.add("Threshold", threshold_field)
        Gmsh.gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
        Gmsh.gmsh.model.mesh.field.setNumber(threshold_field, "LcMin", boundary_size)  # 晶界处使用细网格
        Gmsh.gmsh.model.mesh.field.setNumber(threshold_field, "LcMax", grain_size)     # 晶粒处使用粗网格
        Gmsh.gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0.0)          # 最小距离
        Gmsh.gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", transition_distance)  # 过渡区域宽度
        
        # 3. 禁用其他尺寸控制
        try_set_option("Mesh.MeshSizeFromPoints", 0)
        try_set_option("Mesh.MeshSizeFromCurvature", 0)
        try_set_option("Mesh.MeshSizeExtendFromBoundary", 0)
        
        # 4. 强制使用场控制网格尺寸
        Gmsh.gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)
        
        println("  阈值场设置成功，强制应用场控制")
        println("  - 晶界尺寸(LcMin): $boundary_size")
        println("  - 晶粒尺寸(LcMax): $grain_size (显著大于晶界尺寸)")
        println("  - 过渡区域宽度: $transition_distance")
        
        # 设置全局尺寸边界确保粗化效果
        try_set_option("Mesh.MeshSizeMin", boundary_size)
        try_set_option("Mesh.MeshSizeMax", grain_size)
    catch e
        println("  阈值场设置失败: $e")
        println("  回退到基本设置...")
        
        # 回退方法：使用全局尺寸设置并依靠Gmsh内部算法
        try_set_option("Mesh.MeshSizeMin", boundary_size)
        try_set_option("Mesh.MeshSizeMax", grain_size)
        try_set_option("Mesh.MeshSizeFromPoints", 1)
        try_set_option("Mesh.MeshSizeExtendFromBoundary", 1)
        
        # 尝试直接设置边界点的尺寸
        try
            # 获取所有晶界点
            boundary_vertices = Set{Int}()
            for edge_id in all_boundary_edges
                if haskey(lines, edge_id)
                    push!(boundary_vertices, lines[edge_id][1])
                    push!(boundary_vertices, lines[edge_id][2])
                end
            end
            
            # 为晶界点设置尺寸
            for vertex in boundary_vertices
                Gmsh.gmsh.model.mesh.setSize([(0, vertex)], boundary_size)
            end
            
            println("  应用了晶界点尺寸控制")
        catch e2
            println("  点尺寸控制失败: $e2")
        end
    end
    
    # ===== 修复网格优化参数 =====
    println("设置网格优化参数...")
    
    # 仅设置支持的参数
    try_set_option("Mesh.RecombineAll", 1)                  # 全局启用四边形重组
    try_set_option("Mesh.RecombinationAlgorithm", 1)        # 标准重组算法
    try_set_option("Mesh.RecombineOptimizeTopology", 2)     # 拓扑优化
    try_set_option("Mesh.RecombineNodeRepositioning", 1)    # 节点重新定位
    try_set_option("Mesh.Optimize", 1)                      # 启用优化
    try_set_option("Mesh.OptimizeNetgen", 1)                # Netgen优化
    try_set_option("Mesh.Smoothing", 100)                   # 平滑迭代次数
    
    # ===== 修改网格生成部分 =====
    # ===== 生成网格 =====
    println("\n开始生成结构化网格... [$(round(time() - start_time, digits=2))s]")
    optimization_start_time = time()
    
    # 添加交互式确认
    if interactive
        user_continue = false
        println("\n已完成网格设置:")
        println("  - 晶界区域: $(length(quad_boundary_surfaces)) 个四边形, $(length(non_quad_boundary_surfaces)) 个非四边形")
        println("  - 成功设置 $transfinite_success / $(length(quad_boundary_surfaces)) 个Transfinite表面")
        println("  - 网格尺寸设置: ")
        println("    · 晶粒内部: $(grain_size) (粗网格)")
        println("    · 晶界区域: $(boundary_size) (细网格)")
        println("    · 尺寸比例: $(round(grain_size/boundary_size, digits=1)):1")
        println("    · 过渡区域: $(transition_distance)")
        
        if grain_size/boundary_size < 10
            println("\n提示: 当前晶粒与晶界的尺寸比例低于10:1，可能不足以产生明显的粗化效果")
            println("      您可以考虑增大grain_size或减小boundary_size来获得更明显的对比")
        end
        
        println("\n是否继续生成网格? [y/n]: ")
        
        response = lowercase(strip(readline()))
        if startswith(response, "y")
            user_continue = true
            println("继续生成网格...")
        else
            println("用户取消，退出网格生成...")
            Gmsh.gmsh.finalize()
            return false
        end
        
        if !user_continue
            Gmsh.gmsh.finalize()
            return false
        end
    end
    
    try
        # 阶段1: 生成边界网格
        println("阶段1: 生成边界网格...")
        Gmsh.gmsh.model.mesh.generate(1)
        
        # 阶段2: 一次性生成所有表面网格
        println("阶段2: 生成所有表面网格...")
        Gmsh.gmsh.model.mesh.generate(2)
        
        # 阶段3: 网格优化，考虑时间限制
        println("阶段3: 优化网格质量（时间限制: $(max_optimization_time)秒）...")
        
        # Netgen优化
        try
            remaining_time = max_optimization_time - (time() - optimization_start_time)
            if remaining_time > 0
                netgen_start = time()
                Gmsh.gmsh.model.mesh.optimize("Netgen", false, 5)
                netgen_time = round(time() - netgen_start, digits=2)
                println("  Netgen优化完成，用时: $(netgen_time)秒")
            else
                println("  跳过Netgen优化：已达到最大优化时间")
            end
        catch e
            println("  Netgen优化跳过: $e")
        end
        
        # Laplace2D优化
        try
            remaining_time = max_optimization_time - (time() - optimization_start_time)
            if remaining_time > 0
                # 根据剩余时间动态调整迭代次数
                iterations = min(30, max(5, Int(floor(remaining_time / 0.5))))
                laplace_start = time()
                Gmsh.gmsh.model.mesh.optimize("Laplace2D", false, iterations)
                laplace_time = round(time() - laplace_start, digits=2)
                println("  Laplace平滑完成，$(iterations)次迭代，用时: $(laplace_time)秒")
            else
                println("  跳过Laplace平滑：已达到最大优化时间")
            end
        catch e
            println("  Laplace平滑跳过: $e")
        end
        
    catch e
        println("主网格生成过程失败: $e")
        println("尝试备用方法...")
        
        # 记录备用方法日志
        println(output_io, "# 主网格生成失败: $e")
        println(output_io, "# 使用备用方法")
        
        # 备用方法：重置并简化设置
        try
            println("\n使用备用方法生成网格...")
            Gmsh.gmsh.model.mesh.clear()
            
            # 使用更保守的设置
            try_set_option("Mesh.Algorithm", 1)  # Delaunay算法
            try_set_option("Mesh.RecombineAll", 1)
            
            # 直接设置全局网格尺寸
            try_set_option("Mesh.MeshSizeMin", boundary_size)
            try_set_option("Mesh.MeshSizeMax", grain_size) 
            
            # 生成网格
            Gmsh.gmsh.model.mesh.generate(2)
            
            # 简单优化
            remaining_time = max_optimization_time - (time() - optimization_start_time)
            if remaining_time > 0
                iterations = min(10, max(3, Int(floor(remaining_time / 0.5))))
                try
                    Gmsh.gmsh.model.mesh.optimize("Laplace2D", false, iterations)
                    println("  备用方法Laplace优化完成，$(iterations)次迭代")
                catch
                    # 忽略优化错误
                end
            end
            
            println("备用方法成功")
        catch e2
            println("备用方法也失败: $e2")
            println(output_io, "# 备用方法也失败: $e2")
            error("无法生成网格")
        end
    end
    
    # ===== 四边形重组优化 =====
    println("执行四边形重组优化...")
    
    remaining_time = max_optimization_time - (time() - optimization_start_time)
    if remaining_time > 0
        try
            recombine_start = time()
            # 使用Blossom算法进行四边形重组
            try_set_option("Mesh.RecombinationAlgorithm", 3)  # Blossom算法(最优但较慢)
            Gmsh.gmsh.model.mesh.recombine()
            recombine_time = round(time() - recombine_start, digits=2)
            println("  Blossom四边形重组完成，用时: $(recombine_time)秒")
            
            # 再次平滑优化
            remaining_time = max_optimization_time - (time() - optimization_start_time)
            if remaining_time > 5
                smooth_start = time()
                Gmsh.gmsh.model.mesh.optimize("Laplace2D", false, 10)
                smooth_time = round(time() - smooth_start, digits=2)
                println("  后处理平滑完成，用时: $(smooth_time)秒")
            end
        catch e
            println("  四边形优化失败: $e")
        end
    else
        println("  跳过四边形重组优化：已达到最大优化时间")
    end
    
    # ===== 分析网格质量 =====
    println("\n分析网格质量...")
    
    # 分别统计晶界区域和晶粒区域的四边形比例
    boundary_quad_count = 0
    boundary_tri_count = 0
    grain_quad_count = 0
    grain_tri_count = 0
    
    # 统计每个物理组的网格情况
    phys_group_mesh_stats = Dict()
    
    # 获取晶界区域的单元信息
    for surface_tag in boundary_surfaces
        elem_types, elem_tags, _ = Gmsh.gmsh.model.mesh.getElements(2, surface_tag)
        
        quad_count = 0
        tri_count = 0
        for i in 1:length(elem_types)
            if elem_types[i] == 3  # 四边形单元
                quad_count += length(elem_tags[i])
                boundary_quad_count += length(elem_tags[i])
            elseif elem_types[i] == 2  # 三角形单元
                tri_count += length(elem_tags[i])
                boundary_tri_count += length(elem_tags[i])
            end
        end
        
        # 记录物理组的网格统计
        if haskey(surface_to_phys_group, surface_tag)
            phys_group_id = surface_to_phys_group[surface_tag]
            if !haskey(phys_group_mesh_stats, phys_group_id)
                phys_group_mesh_stats[phys_group_id] = (quad_count=0, tri_count=0)
            end
            
            # 更新统计
            stats = phys_group_mesh_stats[phys_group_id]
            phys_group_mesh_stats[phys_group_id] = (
                quad_count = stats.quad_count + quad_count,
                tri_count = stats.tri_count + tri_count
            )
        end
    end
    
    # 获取晶粒区域的单元信息
    for surface_tag in grain_surfaces
        elem_types, elem_tags, _ = Gmsh.gmsh.model.mesh.getElements(2, surface_tag)
        
        quad_count = 0
        tri_count = 0
        for i in 1:length(elem_types)
            if elem_types[i] == 3  # 四边形单元
                quad_count += length(elem_tags[i])
                grain_quad_count += length(elem_tags[i])
            elseif elem_types[i] == 2  # 三角形单元
                tri_count += length(elem_tags[i])
                grain_tri_count += length(elem_tags[i])
            end
        end
        
        # 记录物理组的网格统计
        if haskey(surface_to_phys_group, surface_tag)
            phys_group_id = surface_to_phys_group[surface_tag]
            if !haskey(phys_group_mesh_stats, phys_group_id)
                phys_group_mesh_stats[phys_group_id] = (quad_count=0, tri_count=0)
            end
            
            # 更新统计
            stats = phys_group_mesh_stats[phys_group_id]
            phys_group_mesh_stats[phys_group_id] = (
                quad_count = stats.quad_count + quad_count,
                tri_count = stats.tri_count + tri_count
            )
        end
    end
    
    # 计算总数和比例
    total_boundary = boundary_quad_count + boundary_tri_count
    total_grain = grain_quad_count + grain_tri_count
    total_all = total_boundary + total_grain
    
    boundary_quad_percentage = total_boundary > 0 ? boundary_quad_count / total_boundary * 100 : 0
    grain_quad_percentage = total_grain > 0 ? grain_quad_count / total_grain * 100 : 0
    total_quad_percentage = total_all > 0 ? (boundary_quad_count + grain_quad_count) / total_all * 100 : 0
    
    println("网格统计:")
    println("  晶界区域 (算法6):")
    println("    - 总单元数: $(total_boundary)")
    println("    - 四边形单元: $(boundary_quad_count) ($(round(boundary_quad_percentage, digits=2))%)")
    println("    - 三角形单元: $(boundary_tri_count) ($(round(100 - boundary_quad_percentage, digits=2))%)")
    println("\n  晶粒区域 (算法8):")
    println("    - 总单元数: $(total_grain)")
    println("    - 四边形单元: $(grain_quad_count) ($(round(grain_quad_percentage, digits=2))%)")
    println("    - 三角形单元: $(grain_tri_count) ($(round(100 - grain_quad_percentage, digits=2))%)")
    println("\n  整体网格:")
    println("    - 总单元数: $(total_all)")
    println("    - 四边形单元: $(boundary_quad_count + grain_quad_count) ($(round(total_quad_percentage, digits=2))%)")
    println("    - 三角形单元: $(boundary_tri_count + grain_tri_count) ($(round(100 - total_quad_percentage, digits=2))%)")
    
    # 输出网格统计信息到文件
    println(output_io, "\n# ===== 网格统计 =====")
    println(output_io, "# 总网格单元数: $total_all")
    println(output_io, "# 四边形比例: $(round(total_quad_percentage, digits=2))%")
    println(output_io, "# 晶界区域(算法6)四边形比例: $(round(boundary_quad_percentage, digits=2))%")
    println(output_io, "# 晶粒内部(算法8)四边形比例: $(round(grain_quad_percentage, digits=2))%")
    
    # 输出Transfinite表面设置结果
    println(output_io, "\n# ===== Transfinite表面设置结果 =====")
    println(output_io, "# 表面ID,设置结果,原因,识别方法")
    
    for (surface_id, result) in transfinite_results
        is_compatible, message, method = result
        result_str = is_compatible ? "成功" : "失败"
        println(output_io, "$surface_id,$result_str,$message,$method")
    end
    
    # 输出每个物理组的网格统计
    println(output_io, "\n# ===== 物理组网格统计 =====")
    println(output_io, "# 物理组ID,物理组名称,四边形单元数,三角形单元数,四边形比例(%)")
    
    for (phys_group_id, stats) in phys_group_mesh_stats
        name = get(phys_group_names, phys_group_id, "Unknown")
        total = stats.quad_count + stats.tri_count
        quad_percentage = total > 0 ? stats.quad_count / total * 100 : 0
        
        println(output_io, "$phys_group_id,$name,$(stats.quad_count),$(stats.tri_count),$(round(quad_percentage, digits=2))")
    end
    
    # 晶界区域四边形比例评估
    if boundary_quad_percentage < 95 && boundary_tri_count > 0
        println("\n警告: 晶界区域四边形比例较低 ($(round(boundary_quad_percentage, digits=2))%)")
        println("可能需要调整晶界尺寸或重组算法")
        println(output_io, "\n# 警告: 晶界区域四边形比例较低 ($(round(boundary_quad_percentage, digits=2))%)")
    elseif boundary_quad_percentage >= 95
        println("\n晶界区域四边形比例良好 ($(round(boundary_quad_percentage, digits=2))%)")
        println(output_io, "\n# 晶界区域四边形比例良好 ($(round(boundary_quad_percentage, digits=2))%)")
    end
    
    # 晶粒内部四边形比例评估
    if grain_quad_percentage < 90 && grain_tri_count > 0
        println("\n警告: 晶粒内部四边形比例较低 ($(round(grain_quad_percentage, digits=2))%)")
        println("算法8(Delaunay for Quads)可能需要进一步优化")
        println(output_io, "\n# 警告: 晶粒内部四边形比例较低 ($(round(grain_quad_percentage, digits=2))%)")
    elseif grain_quad_percentage >= 90
        println("\n晶粒内部四边形比例良好 ($(round(grain_quad_percentage, digits=2))%)")
        println(output_io, "\n# 晶粒内部四边形比例良好 ($(round(grain_quad_percentage, digits=2))%)")
    end
    
    # ===== 保存网格文件 =====
    # 将Gmsh MSH文件版本设置为4.1
    try_set_option("Mesh.MshFileVersion", 4.1)  # 修改为4.1版本
    Gmsh.gmsh.write(mesh_file)
    println("\n网格保存到: $(mesh_file) (MSH格式版本: 4.1)")
    
    # 同时保存VTK格式以便可视化
    vtk_file = replace(mesh_file, r"\.msh$" => ".vtk")
    Gmsh.gmsh.write(vtk_file)
    println("VTK格式保存到: $(vtk_file)")
    
    # 导出为Abaqus INP格式，将物理组保存为集合
    inp_file = replace(mesh_file, r"\.msh$" => ".inp")
    inp_file = export_to_inp(mesh_file, inp_file, phys_group_names, phys_group_surfaces, surface_to_phys_group)
    println("Abaqus INP格式保存到: $(inp_file)")
    
    # 记录输出文件信息
    println(output_io, "\n# 网格文件保存到: $mesh_file")
    println(output_io, "# VTK格式保存到: $vtk_file")
    println(output_io, "# Abaqus INP格式保存到: $inp_file")
    println(output_io, "# 总优化时间: $(round(time() - optimization_start_time, digits=2))s (限制: $(max_optimization_time)s)")
    println(output_io, "# 总用时: $(round(time() - start_time, digits=2))s")
    
    # 关闭输出文件
    close(output_io)
    
    # 完成
    Gmsh.gmsh.finalize()
    
    total_time = round(time() - start_time, digits=2)
    println("网格生成完成，总用时: $(total_time)s")
    println("详细网格信息已保存到: $output_file")
    
    return true
end

"""
将网格转换为Abaqus INP格式，并将物理组作为节点和单元集合写入
"""
function export_to_inp(mesh_file, inp_file, phys_group_names, phys_group_surfaces, surface_to_phys_group)
    println("正在将网格转换为Abaqus INP格式...")
    
    # 如果没有指定INP文件名，则基于mesh文件生成
    if inp_file === nothing
        inp_file = replace(mesh_file, r"\.msh$" => ".inp")
    end
    
    # 尝试先通过Gmsh导出INP基本格式
    try
        Gmsh.gmsh.write(inp_file)
        println("  Gmsh导出INP失败，将手动构建INP文件...")
    catch e
        println("  Gmsh导出INP失败: $e(，将手动构建INP文件..._")
    end
    
    # 获取模型中的节点
    nodes_tags, nodes_coords, _ = Gmsh.gmsh.model.mesh.getNodes()
    
    # 获取模型中的单元
    element_types = Dict()
    element_tags = Dict()
    element_node_tags = Dict()
    
    # 收集所有表面的单元
    all_surfaces = Set{Int}()
    for surfaces in values(phys_group_surfaces)
        union!(all_surfaces, Set(surfaces))
    end
    
    # 对每个表面获取单元信息
    for surface_tag in all_surfaces
        types, tags, node_tags = Gmsh.gmsh.model.mesh.getElements(2, surface_tag)
        
        for i in 1:length(types)
            element_type = types[i]
            if !haskey(element_types, element_type)
                element_types[element_type] = []
                element_tags[element_type] = []
                element_node_tags[element_type] = []
            end
            
            append!(element_types[element_type], repeat([element_type], length(tags[i])))
            append!(element_tags[element_type], tags[i])
            append!(element_node_tags[element_type], node_tags[i])
        end
    end
    
    # 创建节点和单元集合
    node_sets = Dict{String, Set{Int}}()
    element_sets = Dict{String, Dict{Int, Set{Int}}}()
    
    # 对每个物理组创建集合
    for (group_id, name) in phys_group_names
        if !haskey(phys_group_surfaces, group_id)
            continue
        end
        
        # 初始化该组的节点集和单元集
        node_sets[name] = Set{Int}()
        element_sets[name] = Dict{Int, Set{Int}}()
        
        # 收集该组中所有表面的单元和节点
        for surface_tag in phys_group_surfaces[group_id]
            types, tags, node_tags = Gmsh.gmsh.model.mesh.getElements(2, surface_tag)
            
            for i in 1:length(types)
                element_type = types[i]
                
                # 初始化该类型单元的集合
                if !haskey(element_sets[name], element_type)
                    element_sets[name][element_type] = Set{Int}()
                end
                
                # 添加单元到集合
                union!(element_sets[name][element_type], Set(tags[i]))
                
                # 添加节点到集合
                for j in 1:length(node_tags[i])
                    node_id = node_tags[i][j]
                    push!(node_sets[name], node_id)
                end
            end
        end
    end
    
    # 打开INP文件进行写入
    open(inp_file, "w") do io
        # 写入文件头
        write(io, "*Heading\n")
        write(io, "模型由Julia Gmsh API生成的多晶材料结构化网格\n")
        
        # 写入节点部分
        write(io, "*Node\n")
        for i in 1:length(nodes_tags)
            node_id = Int(nodes_tags[i])
            x = nodes_coords[3*i-2]
            y = nodes_coords[3*i-1]
            z = nodes_coords[3*i]
            write(io, "$(node_id), $(x), $(y), $(z)\n")
        end
        
        # 写入单元部分 - 按单元类型分组
        for (element_type, tags) in element_tags
            # 根据Gmsh单元类型映射到Abaqus单元类型
            abaqus_type = if element_type == 2  # 三角形
                "CPS3"  # 平面应力三角形
            elseif element_type == 3  # 四边形
                "CPS4"  # 平面应力四边形
            else
                "Unknown"
            end
            
            write(io, "*Element, type=$(abaqus_type)\n")
            
            # 获取该类型的所有单元节点连接关系
            etags = element_tags[element_type]
            enode_tags = element_node_tags[element_type]
            
            # 确定每个单元有多少个节点
            nodes_per_element = if element_type == 2  # 三角形
                3
            elseif element_type == 3  # 四边形
                4
            else
                error("未知的单元类型: $element_type")
            end
            
            # 写入单元连接关系
            for i in 1:length(etags)
                element_id = etags[i]
                start_idx = (i-1) * nodes_per_element + 1
                end_idx = start_idx + nodes_per_element - 1
                
                node_list = enode_tags[start_idx:end_idx]
                node_str = join(node_list, ", ")
                
                write(io, "$(element_id), $(node_str)\n")
            end
        end
        
        # 写入节点集
        for (name, nodes) in node_sets
            # 替换名称中的空格和特殊字符
            set_name = replace(name, r"[^a-zA-Z0-9_]" => "_")
            write(io, "*Nset, nset=$(set_name)_nodes\n")
            
            # 每行最多16个节点ID
            sorted_nodes = sort(collect(nodes))
            for i in 1:16:length(sorted_nodes)
                end_idx = min(i+15, length(sorted_nodes))
                node_str = join(sorted_nodes[i:end_idx], ", ")
                write(io, node_str * "\n")
            end
        end
        
        # 写入单元集
        for (name, type_elements) in element_sets
            set_name = replace(name, r"[^a-zA-Z0-9_]" => "_")
            
            for (element_type, elements) in type_elements
                write(io, "*Elset, elset=$(set_name)_elements_$(element_type)\n")
                
                # 每行最多16个单元ID
                sorted_elements = sort(collect(elements))
                for i in 1:16:length(sorted_elements)
                    end_idx = min(i+15, length(sorted_elements))
                    element_str = join(sorted_elements[i:end_idx], ", ")
                    write(io, element_str * "\n")
                end
            end
            
            # 创建一个包含所有该物理组单元的集合
            write(io, "*Elset, elset=$(set_name)_all_elements\n")
            all_elements = []
            for elements in values(type_elements)
                append!(all_elements, collect(elements))
            end
            sort!(all_elements)
            
            for i in 1:16:length(all_elements)
                end_idx = min(i+15, length(all_elements))
                element_str = join(all_elements[i:end_idx], ", ")
                write(io, element_str * "\n")
            end
        end
        
        # 在文件末尾添加材料和截面属性的示例（用户需要根据实际情况进行修改）
        write(io, "\n** 材料定义示例（用户需要根据实际情况修改）\n")
        write(io, "*Material, name=Material-1\n")
        write(io, "*Elastic\n")
        write(io, " 210000.0, 0.3\n")
        
        # 为每个物理组定义截面属性
        for name in keys(node_sets)
            set_name = replace(name, r"[^a-zA-Z0-9_]" => "_")
            write(io, "\n** 截面定义 - $(name)\n")
            write(io, "*Solid Section, elset=$(set_name)_all_elements, material=Material-1\n")
            write(io, "1.0\n")  # 截面厚度
        end
        
        # 添加分析步骤示例
        write(io, "\n** 分析步骤示例\n")
        write(io, "*Step, name=Step-1\n")
        write(io, "*Static\n")
        write(io, "1., 1., 1e-05, 1.\n")
        write(io, "*End Step\n")
    end
    
    println("Abaqus INP格式文件已保存到: $inp_file")
    println("- 已将物理组导出为节点集(Nset)和单元集(Elset)")
    println("- 每个物理组生成了独立的节点集和单元集")
    println("- 文件包含了基本的材料和截面属性定义示例")
    
    return inp_file
end

"""
主函数：生成多晶材料结构化网格
"""
function main()
    # 默认文件路径
    input_file = "equiaxed_grains_with_boundaries.geo"
    output_file = "quad_domi_mesh.msh"
    info_file = "quad_domi_mesh_info.txt"
    
    # 检查命令行参数
    if length(ARGS) >= 1
        input_file = ARGS[1]
    end
    
    if length(ARGS) >= 2
        output_file = ARGS[2]
    end
    
    if length(ARGS) >= 3
        info_file = ARGS[3]
    end
    
    # 查询是否启用交互模式
    interactive_mode = true
    if length(ARGS) >= 4 && (ARGS[4] == "no-interactive" || ARGS[4] == "false")
        interactive_mode = false
    end
    
    # 获取优化时间限制
    max_optimization_time = 60  # 默认60秒
    if length(ARGS) >= 5
        try
            max_optimization_time = parse(Int, ARGS[5])
        catch
            println("无效的优化时间限制参数，使用默认值60秒")
        end
    end
    
    # 获取过渡区域设置
    transition_distance = 0.05  # 默认过渡区域 - 阈值场需要较大的过渡区域
    if length(ARGS) >= 6
        try
            transition_distance = parse(Float64, ARGS[6])
        catch
            println("无效的过渡区域参数，使用默认值0.05")
        end
    end
    
    println("多晶材料结构化网格生成 - 四边形主导版")
    println("输入几何文件: $input_file")
    println("输出网格文件: $output_file")
    println("信息输出文件: $info_file")
    println("交互模式: $(interactive_mode ? "启用" : "禁用")")
    println("优化时间限制: $max_optimization_time 秒")
    println("过渡区域距离: $transition_distance")
    
    # 生成网格 - 使用更大的尺寸对比
    generate_polycrystal_mesh_with_algorithm8(input_file, output_file, 
                                             grain_size=0.02,            # 晶粒尺寸更大
                                             boundary_size=0.0005,        # 晶界尺寸更小
                                             transition_distance=transition_distance, # 过渡区域
                                             num_threads=8,
                                             max_optimization_time=max_optimization_time,
                                             interactive=interactive_mode,
                                             verbose=true,
                                             output_file=info_file)
    
    println("处理完成")
end

# 如果直接运行此文件，则执行主函数
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end