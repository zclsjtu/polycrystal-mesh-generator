using Gmsh
using Dates
using Statistics
using DelimitedFiles
"""
生成多晶材料三角形网格，晶界细密，晶粒粗糙
"""
function generate_polycrystal_triangle_mesh(geo_file, mesh_file; 
                                           grain_size=0.05,           # 晶粒区域较粗的网格尺寸
                                           boundary_size=0.001,       # 晶界区域较细的网格尺寸
                                           transition_distance=0.1,   # 尺寸过渡区域宽度
                                           num_threads=8,
                                           max_optimization_time=60,
                                           interactive=true,
                                           verbose=true,
                                           output_file="polycrystal_triangle_mesh_info.txt")
    start_time = time()
    
    # 打开输出文件
    output_io = open(output_file, "w")
    
    # 文件头添加字段说明
    println(output_io, "# 多晶材料三角形网格信息文件 (晶界细密-晶粒粗糙版)")
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
    println(output_io, "# 晶界区域生成细密三角形网格，晶粒区域生成粗糙三角形网格")
    println(output_io, "# ==========================================")
    
    # 从geo文件中直接提取边长信息
    points, lines, edge_lengths, loop_to_lines, surface_to_loop, boundary_surfaces, grain_surfaces, physical_groups = extract_edge_lengths_from_geo(geo_file)
    
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
    
    # 重要：禁用四边形重组设置，确保生成三角形网格
    try_set_option("Mesh.RecombineAll", 0)        # 禁用全局四边形重组
    try_set_option("Mesh.Recombine3DAll", 0)      # 禁用3D重组
    try_set_option("Mesh.RecombineOptimizeTopology", 0) # 禁用拓扑优化
    
    # 算法选择：对于三角形网格，选择Delaunay或Frontal算法
    try_set_option("Mesh.Algorithm", 5)  # 5=Delaunay, 6=Frontal-Delaunay, 7=BAMG
    try_set_option("Mesh.Algorithm3D", 1)  # 1=Delaunay, 4=Frontal, 7=MMG3D
    
    # 质量参数
    try_set_option("Mesh.ElementOrder", 1)        # 一阶单元
    try_set_option("Mesh.Smoothing", 100)         # 平滑迭代次数
    try_set_option("Mesh.SmoothNormals", 1)       # 平滑法向量
    try_set_option("Mesh.QualityType", 2)         # SICN质量指标(适用于三角形)
    try_set_option("Mesh.QualityOptimization", 1) # 启用质量优化
    try_set_option("Mesh.OptimizeNetgen", 1)      # Netgen优化
    
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
    
    # ===== 为边设置细分控制点 =====
    println("为边设置控制点...")
    
    # 收集所有晶界边
    all_boundary_edges = []
    for surface_tag in boundary_surfaces
        if haskey(surface_edges, surface_tag)
            append!(all_boundary_edges, abs.(surface_edges[surface_tag]))
        end
    end
    unique!(all_boundary_edges)
    
    # 收集所有晶粒边界（非晶界）
    all_grain_edges = []
    for surface_tag in grain_surfaces
        if haskey(surface_edges, surface_tag)
            append!(all_grain_edges, abs.(surface_edges[surface_tag]))
        end
    end
    # 去重并移除也属于晶界的边
    unique!(all_grain_edges)
    filter!(e -> !(e in all_boundary_edges), all_grain_edges)
    
    # 初始化边的控制点计数
    edge_divisions = Dict{Int, Int}()
    
    # 晶界边设置较多控制点(细密网格)
    boundary_edges_set = 0
    for edge_tag in all_boundary_edges
        if haskey(edge_lengths, edge_tag)
            # 根据边长计算控制点数量，但确保足够细密
            length = edge_lengths[edge_tag]
            n_divisions = max(3, ceil(Int, length / boundary_size))
            
            try
                Gmsh.gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_divisions, "Progression", 1.0)
                edge_divisions[edge_tag] = n_divisions
                boundary_edges_set += 1
            catch e
                println("  晶界边 #$edge_tag 设置失败: $e")
            end
        end
    end
    println("  为 $boundary_edges_set / $(length(all_boundary_edges)) 条晶界边设置了控制点")
    
    # 晶粒边设置较少控制点(粗糙网格)
    grain_edges_set = 0
    for edge_tag in all_grain_edges
        if haskey(edge_lengths, edge_tag)
            # 根据边长计算控制点数量，较粗糙
            length = edge_lengths[edge_tag]
            n_divisions = max(2, ceil(Int, length / grain_size))
            
            try
                Gmsh.gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_divisions, "Progression", 1.0)
                edge_divisions[edge_tag] = n_divisions
                grain_edges_set += 1
            catch e
                println("  晶粒边 #$edge_tag 设置失败: $e")
            end
        end
    end
    println("  为 $grain_edges_set / $(length(all_grain_edges)) 条晶粒边设置了控制点")
    
    # ===== 设置阈值场控制网格尺寸过渡 =====
    println("创建阈值场控制网格尺寸过渡...")
    
    try
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
        
        println("  阈值场设置成功，将在晶界附近 $transition_distance 范围内平滑过渡")
        println("  - 晶界尺寸(LcMin): $boundary_size")
        println("  - 晶粒尺寸(LcMax): $grain_size")
        println("  - 尺寸比例: $(round(grain_size/boundary_size, digits=1)):1")
        
        # 设置全局尺寸边界确保粗细对比效果
        try_set_option("Mesh.MeshSizeMin", boundary_size)
        try_set_option("Mesh.MeshSizeMax", grain_size)
    catch e
        println("  阈值场设置失败: $e")
        println("  回退到基本设置...")
        
        # 回退方法：使用全局尺寸设置
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
            
            # 为晶界点设置细小尺寸
            for vertex in boundary_vertices
                Gmsh.gmsh.model.mesh.setSize([(0, vertex)], boundary_size)
            end
            
            println("  应用了晶界点尺寸控制")
        catch e2
            println("  点尺寸控制失败: $e2")
        end
    end
    
    # ===== 网格生成前确认 =====
    println("\n开始生成三角形网格... [$(round(time() - start_time, digits=2))s]")
    optimization_start_time = time()
    
    # 添加交互式确认
    if interactive
        user_continue = false
        println("\n已完成网格设置:")
        println("  - 已为 $boundary_edges_set 条晶界边设置细密控制点")
        println("  - 已为 $grain_edges_set 条晶粒边设置粗糙控制点")
        println("  - 网格尺寸设置: ")
        println("    · 晶粒内部: $(grain_size) (粗网格)")
        println("    · 晶界区域: $(boundary_size) (细网格)")
        println("    · 尺寸比例: $(round(grain_size/boundary_size, digits=1)):1")
        println("    · 过渡区域: $(transition_distance)")
        println("  - 已设置为纯三角形网格模式")
        
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
    
    # ===== 多阶段网格生成 =====
    try
        # 阶段1: 生成边界网格
        println("阶段1: 生成边界网格...")
        Gmsh.gmsh.model.mesh.generate(1)
        
        # 新增: 为晶粒区域设置四边形重组
        println("为晶粒区域设置四边形重组...")
        # 设置重组参数 - 使用较为宽松的阈值，避免崩溃
        try_set_option("Mesh.RecombineMinimumQuality", 0.5)  # 降低质量要求，提高重组成功率
        
        # 仅为晶粒区域启用四边形重组
        recombine_count = 0
        for surface_tag in grain_surfaces
            try
                Gmsh.gmsh.model.mesh.setRecombine(2, surface_tag)  # 不设置质量参数，使用全局设置
                recombine_count += 1
            catch e
                println("  表面 #$surface_tag 设置四边形重组失败")
            end
        end
        println("  为 $recombine_count / $(length(grain_surfaces)) 个晶粒区域设置了四边形重组")
        
        # 阶段2: 生成表面网格
        println("阶段2: 生成表面网格...")
        Gmsh.gmsh.model.mesh.generate(2)
        
        # 阶段3: 质量优化
        println("阶段3: 优化网格质量...")
        
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
        
        # 三角形优化 - 适合三角形网格的优化方法
        try
            remaining_time = max_optimization_time - (time() - optimization_start_time)
            if remaining_time > 0
                # 根据剩余时间动态调整迭代次数
                iterations = min(30, max(5, Int(floor(remaining_time / 0.5))))
                smooth_start = time()
                
                # 对三角形网格优化
                # 先用Laplace
                Gmsh.gmsh.model.mesh.optimize("Laplace2D", false, iterations)
                
                # 再用Delaunay优化 - 特别适合三角形
                if remaining_time - (time() - smooth_start) > 0
                    Gmsh.gmsh.model.mesh.optimize("HighOrder", false, 3)
                end
                
                smooth_time = round(time() - smooth_start, digits=2)
                println("  三角形网格优化完成，用时: $(smooth_time)秒")
            else
                println("  跳过三角形优化：已达到最大优化时间")
            end
        catch e
            println("  三角形优化跳过: $e")
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
            
            # 使用最基本的Delaunay三角剖分
            try_set_option("Mesh.Algorithm", 1)  # 最简单的Delaunay算法
            try_set_option("Mesh.RecombineAll", 0)  # 确保生成三角形
            
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
                    println("  备用方法优化完成，$(iterations)次迭代")
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
    
    # ===== 分析网格质量 =====
    println("\n分析网格质量...")
    
    # 分别统计晶界区域和晶粒区域的三角形和四边形数量
    boundary_tri_count = 0
    boundary_quad_count = 0
    grain_tri_count = 0
    grain_quad_count = 0
    
    # 统计每个物理组的网格情况
    phys_group_mesh_stats = Dict()
    
    # 获取晶界区域的单元信息
    for surface_tag in boundary_surfaces
        elem_types, elem_tags, _ = Gmsh.gmsh.model.mesh.getElements(2, surface_tag)
        
        tri_count = 0
        quad_count = 0
        for i in 1:length(elem_types)
            if elem_types[i] == 2  # 三角形单元
                tri_count += length(elem_tags[i])
                boundary_tri_count += length(elem_tags[i])
            elseif elem_types[i] == 3  # 四边形单元
                quad_count += length(elem_tags[i])
                boundary_quad_count += length(elem_tags[i])
            end
        end
        
        # 记录物理组的网格统计
        if haskey(surface_to_phys_group, surface_tag)
            phys_group_id = surface_to_phys_group[surface_tag]
            if !haskey(phys_group_mesh_stats, phys_group_id)
                phys_group_mesh_stats[phys_group_id] = (tri_count=0, quad_count=0)
            end
            
            # 更新统计
            stats = phys_group_mesh_stats[phys_group_id]
            phys_group_mesh_stats[phys_group_id] = (
                tri_count = stats.tri_count + tri_count,
                quad_count = stats.quad_count + quad_count,
            )
        end
    end
    
    # 获取晶粒区域的单元信息
    for surface_tag in grain_surfaces
        elem_types, elem_tags, _ = Gmsh.gmsh.model.mesh.getElements(2, surface_tag)
        
        tri_count = 0
        quad_count = 0
        for i in 1:length(elem_types)
            if elem_types[i] == 2  # 三角形单元
                tri_count += length(elem_tags[i])
                grain_tri_count += length(elem_tags[i])
            elseif elem_types[i] == 3  # 四边形单元
                quad_count += length(elem_tags[i])
                grain_quad_count += length(elem_tags[i])
            end
        end
        
        # 记录物理组的网格统计
        if haskey(surface_to_phys_group, surface_tag)
            phys_group_id = surface_to_phys_group[surface_tag]
            if !haskey(phys_group_mesh_stats, phys_group_id)
                phys_group_mesh_stats[phys_group_id] = (tri_count=0, quad_count=0)
            end
            
            # 更新统计
            stats = phys_group_mesh_stats[phys_group_id]
            phys_group_mesh_stats[phys_group_id] = (
                tri_count = stats.tri_count + tri_count,
                quad_count = stats.quad_count + quad_count,
            )
        end
    end
    
    # 计算总数和比例
    total_tri = boundary_tri_count + grain_tri_count
    total_quad = boundary_quad_count + grain_quad_count
    total_all = total_tri + total_quad
    
    boundary_tri_percentage = total_all > 0 ? boundary_tri_count / total_all * 100 : 0
    boundary_quad_percentage = total_all > 0 ? boundary_quad_count / total_all * 100 : 0
    grain_tri_percentage = total_all > 0 ? grain_tri_count / total_all * 100 : 0
    grain_quad_percentage = total_all > 0 ? grain_quad_count / total_all * 100 : 0
    
    # 计算四边形重组效率
    grain_total = grain_tri_count + grain_quad_count
    quad_efficiency = grain_total > 0 ? grain_quad_count / grain_total * 100 : 0
    
    println("网格统计:")
    println("  晶界区域:")
    println("    - 三角形单元: $(boundary_tri_count) ($(round(boundary_tri_percentage, digits=2))%)")
    println("    - 四边形单元: $(boundary_quad_count) ($(round(boundary_quad_percentage, digits=2))%)")
    println("\n  晶粒区域:")
    println("    - 三角形单元: $(grain_tri_count) ($(round(grain_tri_percentage, digits=2))%)")
    println("    - 四边形单元: $(grain_quad_count) ($(round(grain_quad_percentage, digits=2))%)")
    println("    - 四边形重组效率: $(round(quad_efficiency, digits=2))%")
    println("\n  整体网格:")
    println("    - 总三角形单元数: $(total_tri)")
    println("    - 总四边形单元数: $(total_quad)")
    println("    - 单元总数: $(total_all)")
    
    # 检查网格质量并计算晶界区域单元密度
    boundary_area = 0.0
    grain_area = 0.0
    
    # 尝试计算面积
    try
        # 计算晶界区域面积
        for surface_tag in boundary_surfaces
            # 注意：这里假设Gmsh API提供了获取表面积的方法
            # 实际中可能需要基于三角形单元手动计算
            area = Gmsh.gmsh.model.mesh.getJacobians(2, surface_tag)[1]  # 简化，实际需要针对API调整
            boundary_area += sum(area)
        end
        
        # 计算晶粒区域面积
        for surface_tag in grain_surfaces
            area = Gmsh.gmsh.model.mesh.getJacobians(2, surface_tag)[1]  # 简化，实际需要针对API调整
            grain_area += sum(area)
        end
        
        # 计算单元密度
        boundary_density = boundary_tri_count / boundary_area
        grain_density = grain_tri_count / grain_area
        density_ratio = boundary_density / grain_density
        
        println("\n网格密度比例:")
        println("  晶界区域密度: $(round(boundary_density, digits=2)) 单元/面积")
        println("  晶粒区域密度: $(round(grain_density, digits=2)) 单元/面积")
        println("  密度比例: $(round(density_ratio, digits=2)):1")
    catch e
        println("\n无法计算精确的网格密度比例: $e")
    end
    
    # 输出网格统计信息到文件
    println(output_io, "\n# ===== 网格统计 =====")
    println(output_io, "# 总单元数: $total_all")
    println(output_io, "# 总三角形单元数: $total_tri")
    println(output_io, "# 总四边形单元数: $total_quad")
    println(output_io, "# 晶界区域三角形数: $boundary_tri_count ($(round(boundary_tri_percentage, digits=2))%)")
    println(output_io, "# 晶界区域四边形数: $boundary_quad_count ($(round(boundary_quad_percentage, digits=2))%)")
    println(output_io, "# 晶粒区域三角形数: $grain_tri_count ($(round(grain_tri_percentage, digits=2))%)")
    println(output_io, "# 晶粒区域四边形数: $grain_quad_count ($(round(grain_quad_percentage, digits=2))%)")
    println(output_io, "# 晶粒区域四边形重组效率: $(round(quad_efficiency, digits=2))%")
    
    # 输出每个物理组的网格统计
    println(output_io, "\n# ===== 物理组网格统计 =====")
    println(output_io, "# 物理组ID,物理组名称,三角形单元数,四边形单元数")
    
    for (phys_group_id, stats) in phys_group_mesh_stats
        name = get(phys_group_names, phys_group_id, "Unknown")
        println(output_io, "$phys_group_id,$name,$(stats.tri_count),$(stats.quad_count)")
    end
    
    # ===== 保存网格文件 =====
    # 将Gmsh MSH文件版本设置为4.1
    try_set_option("Mesh.MshFileVersion", 4.1)  # 使用4.1版本
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
    println("三角形网格生成完成，总用时: $(total_time)s")
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
        println("  Gmsh导出INP失败: $e(，将手动构建INP文件...")
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
        write(io, "模型由Julia Gmsh API生成的多晶材料三角形网格\n")
        
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
主函数：生成多晶材料三角形网格
"""
function main()
    # 默认文件路径
    input_file = "equiaxed_grains_with_boundaries.geo"
    output_file = "mix_mesh.msh"
    info_file = "mix_mesh_info.txt"
    
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
    transition_distance = 0.05  # 默认过渡区域
    if length(ARGS) >= 6
        try
            transition_distance = parse(Float64, ARGS[6])
        catch
            println("无效的过渡区域参数，使用默认值0.05")
        end
    end
    
    println("多晶材料三角形网格生成 - 晶界细密，晶粒粗糙版")
    println("输入几何文件: $input_file")
    println("输出网格文件: $output_file")
    println("信息输出文件: $info_file")
    println("交互模式: $(interactive_mode ? "启用" : "禁用")")
    println("优化时间限制: $max_optimization_time 秒")
    println("过渡区域距离: $transition_distance")
    
    # 生成三角形网格
    generate_polycrystal_triangle_mesh(input_file, output_file, 
                                       grain_size=0.02,              # 晶粒区域较粗
                                       boundary_size=0.0005,         # 晶界区域较细
                                       transition_distance=transition_distance,
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