using LinearAlgebra
using Plots
using DelimitedFiles

"""
解析geo文件提取顶点和晶粒
"""
function read_geo_file(filename)
    # 读取文件内容
    if !isfile(filename)
        error("文件不存在: $filename")
    end
    
    points = Dict{Int, Vector{Float64}}()  # 点ID -> 坐标
    lines = Dict{Int, Vector{Int}}()       # 线ID -> [起点ID, 终点ID]
    line_loops = Dict{Int, Vector{Int}}()  # 线环ID -> [线ID数组]
    surfaces = Dict{Int, Int}()            # 表面ID -> 线环ID
    
    file_content = read(filename, String)
    
    # 提取点
    for m in eachmatch(r"Point\s*\(\s*(\d+)\s*\)\s*=\s*\{\s*([^,]+)\s*,\s*([^,]+)", file_content)
        id = parse(Int, m.captures[1])
        x = parse(Float64, m.captures[2])
        y = parse(Float64, m.captures[3])
        points[id] = [x, y]
    end
    
    # 提取线
    for m in eachmatch(r"Line\s*\(\s*(\d+)\s*\)\s*=\s*\{\s*(\d+)\s*,\s*(\d+)\s*\}", file_content)
        id = parse(Int, m.captures[1])
        start_id = parse(Int, m.captures[2])
        end_id = parse(Int, m.captures[3])
        lines[id] = [start_id, end_id]
    end
    
    # 提取线环
    for m in eachmatch(r"Line\s*Loop\s*\(\s*(\d+)\s*\)\s*=\s*\{([^}]+)\}", file_content)
        id = parse(Int, m.captures[1])
        line_ids_str = split(m.captures[2], ",")
        line_ids = [parse(Int, strip(s)) for s in line_ids_str]
        line_loops[id] = line_ids
    end
    
    # 提取表面
    for m in eachmatch(r"(Plane\s*Surface|Surface)\s*\(\s*(\d+)\s*\)\s*=\s*\{([^}]+)\}", file_content)
        id = parse(Int, m.captures[2])
        loop_ids_str = split(m.captures[3], ",")
        if !isempty(loop_ids_str)
            loop_id = parse(Int, strip(loop_ids_str[1]))
            surfaces[id] = loop_id
        end
    end
    
    println("解析文件: $filename")
    println("  找到 $(length(points)) 个点")
    println("  找到 $(length(lines)) 条线")
    println("  找到 $(length(line_loops)) 个线环")
    println("  找到 $(length(surfaces)) 个表面")
    
    return points, lines, line_loops, surfaces
end

"""
获取从线环提取的有序顶点列表和顶点ID
"""
function get_ordered_vertices_from_loop(loop_id, line_loops, lines, points)
    if !haskey(line_loops, loop_id)
        return Int[], Vector{Float64}[]
    end
    
    line_ids = line_loops[loop_id]
    if isempty(line_ids)
        return Int[], Vector{Float64}[]
    end
    
    # 创建连接的线段列表
    segments = []
    for line_id in line_ids
        abs_id = abs(line_id)
        if haskey(lines, abs_id)
            segment = copy(lines[abs_id])
            if line_id < 0  # 负号表示反向
                segment = reverse(segment)
            end
            push!(segments, segment)
        end
    end
    
    if isempty(segments)
        return Int[], Vector{Float64}[]
    end
    
    # 创建有序的点列表
    ordered_segments = []
    remaining_segments = copy(segments)
    
    # 从第一个线段开始
    push!(ordered_segments, popfirst!(remaining_segments))
    
    # 继续添加连接的线段
    while !isempty(remaining_segments)
        current_end = ordered_segments[end][2]
        found = false
        
        for (i, segment) in enumerate(remaining_segments)
            if segment[1] == current_end
                push!(ordered_segments, segment)
                deleteat!(remaining_segments, i)
                found = true
                break
            elseif segment[2] == current_end
                push!(ordered_segments, reverse(segment))
                deleteat!(remaining_segments, i)
                found = true
                break
            end
        end
        
        if !found
            # 如果没有找到连接的线段，试着从起点连接
            current_start = ordered_segments[1][1]
            
            for (i, segment) in enumerate(remaining_segments)
                if segment[2] == current_start
                    prepend!(ordered_segments, [segment])
                    deleteat!(remaining_segments, i)
                    found = true
                    break
                elseif segment[1] == current_start
                    prepend!(ordered_segments, [reverse(segment)])
                    deleteat!(remaining_segments, i)
                    found = true
                    break
                end
            end
            
            if !found && !isempty(remaining_segments)
                # 如果仍然找不到连接的线段，假设有断开的多边形
                # 简单地添加下一个未使用的线段并继续
                push!(ordered_segments, popfirst!(remaining_segments))
            end
        end
    end
    
    # 从有序线段中提取顶点ID
    ordered_vertex_ids = Int[]
    for segment in ordered_segments
        push!(ordered_vertex_ids, segment[1])
    end
    
    # 添加最后一个线段的终点
    if !isempty(ordered_segments)
        push!(ordered_vertex_ids, ordered_segments[end][2])
    end
    
    # 如果是闭合多边形，移除重复的尾部顶点
    if length(ordered_vertex_ids) > 1 && ordered_vertex_ids[1] == ordered_vertex_ids[end]
        pop!(ordered_vertex_ids)
    end
    
    # 获取顶点坐标
    ordered_vertices = Vector{Float64}[]
    for id in ordered_vertex_ids
        if haskey(points, id)
            push!(ordered_vertices, points[id])
        end
    end
    
    return ordered_vertex_ids, ordered_vertices
end

"""
计算多边形质心
"""
function calculate_centroid(vertices)
    if isempty(vertices)
        return [0.0, 0.0]
    end
    
    x_sum = sum(v[1] for v in vertices)
    y_sum = sum(v[2] for v in vertices)
    
    return [x_sum/length(vertices), y_sum/length(vertices)]
end

"""
计算边的法向量（指向多边形内部）
"""
function calculate_edge_normal(p1, p2, centroid)
    edge = p2 - p1
    # 垂直于边的向量
    normal = [edge[2], -edge[1]]
    normal_length = norm(normal)
    
    if normal_length > 0
        normal = normal / normal_length
    else
        return [0.0, 0.0]
    end
    
    # 确保法向量指向多边形内部
    mid_point = (p1 + p2) / 2
    to_center = centroid - mid_point
    
    if dot(normal, to_center) < 0
        normal = -normal
    end
    
    return normal
end

"""
计算两条线的交点
"""
function calculate_intersection(p1, v1, p2, v2)
    # 计算 p1 + t1*v1 = p2 + t2*v2 的交点
    # 重组为 t1*v1 - t2*v2 = p2 - p1
    # 写成矩阵形式：[v1 -v2] * [t1; t2] = p2 - p1
    
    A = [v1[1] -v2[1]; v1[2] -v2[2]]
    b = [p2[1] - p1[1]; p2[2] - p1[2]]
    
    # 检查行列式是否接近零（平行线）
    det_A = A[1,1]*A[2,2] - A[1,2]*A[2,1]
    if abs(det_A) < 1e-10
        return nothing
    end
    
    # 解方程组求t1
    t1 = (A[2,2]*b[1] - A[1,2]*b[2]) / det_A
    
    # 计算交点
    intersection = p1 + t1 * v1
    
    return intersection
end

"""
使用边平移算法生成偏移晶粒
边平移算法步骤：
1. 沿着每条边的法向量方向平移边
2. 计算相邻平移边的交点，这些交点就是偏移多边形的顶点
3. 确保每个老顶点和对应的新顶点正确匹配
"""
function generate_offset_grain_edge_shifting(vertices, offset_distance)
    n = length(vertices)
    if n < 3
        return Vector{Float64}[], Dict{Int, Int}()
    end
    
    centroid = calculate_centroid(vertices)
    
    # 计算每条边的偏移线
    offset_lines = []  # [(起点, 方向向量), ...]
    
    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[i % n + 1]
        
        # 计算边的法向量（指向内部）
        normal = calculate_edge_normal(p1, p2, centroid)
        
        # 偏移边的两个端点
        offset_p1 = p1 + offset_distance * normal
        offset_p2 = p2 + offset_distance * normal
        
        # 偏移边的方向向量
        direction = offset_p2 - offset_p1
        
        push!(offset_lines, (offset_p1, direction))
    end
    
    # 计算相邻偏移线的交点，形成偏移多边形的顶点
    offset_vertices = Vector{Float64}[]
    
    # 顶点对应关系映射：原顶点索引 -> 新顶点索引
    vertex_mapping = Dict{Int, Int}()
    
    for i in 1:n
        current_line = offset_lines[i]
        next_line = offset_lines[i % n + 1]
        
        # 计算交点
        intersection = calculate_intersection(
            current_line[1], current_line[2],
            next_line[1], next_line[2]
        )
        
        if intersection !== nothing
            push!(offset_vertices, intersection)
            # 记录当前原始顶点(i%n+1)与新顶点(length(offset_vertices))的对应关系
            vertex_mapping[i % n + 1] = length(offset_vertices)
        else
            # 如果没有交点（可能是平行线），使用当前线的端点
            fallback_point = current_line[1] + 0.5 * current_line[2]
            push!(offset_vertices, fallback_point)
            vertex_mapping[i % n + 1] = length(offset_vertices)
        end
    end
    
    return offset_vertices, vertex_mapping
end

"""
提取所有晶粒并生成偏移晶粒
"""
function extract_and_offset_grains(points, lines, line_loops, surfaces, offset_distance)
    grains = Dict{Int, Vector{Vector{Float64}}}()
    grain_vertex_ids = Dict{Int, Vector{Int}}()
    offset_grains = Dict{Int, Vector{Vector{Float64}}}()
    vertex_mappings = Dict{Int, Dict{Int, Int}}()
    
    for (surface_id, loop_id) in surfaces
        vertex_ids, vertices = get_ordered_vertices_from_loop(loop_id, line_loops, lines, points)
        
        if length(vertices) >= 3
            grains[surface_id] = vertices
            grain_vertex_ids[surface_id] = vertex_ids
            
            # 使用边平移算法生成偏移晶粒
            offset_vertices, vertex_mapping = generate_offset_grain_edge_shifting(vertices, offset_distance)
            
            if length(offset_vertices) >= 3
                offset_grains[surface_id] = offset_vertices
                vertex_mappings[surface_id] = vertex_mapping
            end
        end
    end
    
    return grains, grain_vertex_ids, offset_grains, vertex_mappings
end

"""
生成晶界模型
"""
function generate_grain_boundary_model(input_file, output_file, offset_distance=0.01)
    # 解析geo文件
    points, lines, line_loops, surfaces = read_geo_file(input_file)
    
    # 提取晶粒并生成偏移晶粒
    grains, grain_vertex_ids, offset_grains, vertex_mappings = 
        extract_and_offset_grains(points, lines, line_loops, surfaces, offset_distance)
    
    println("提取了 $(length(grains)) 个晶粒")
    println("生成了 $(length(offset_grains)) 个偏移晶粒")
    
    # 写入晶界模型
    write_grain_boundary_model(output_file, grains, grain_vertex_ids, offset_grains, vertex_mappings, points)
    
    # 可视化
    visualize_grain_boundaries(grains, offset_grains, vertex_mappings, "grain_boundaries.png")
    
    return grains, offset_grains
end

"""
写入晶界模型
"""
function write_grain_boundary_model(filename, grains, grain_vertex_ids, offset_grains, vertex_mappings, original_points)
    open(filename, "w") do file
        write(file, "// 晶界模型 - 使用边平移算法生成偏移顶点并确保正确的顶点对应关系\n\n")
        write(file, "lc = 0.05;\n\n")
        
        # ============ 写入点 ============
        write(file, "// 原始点\n")
        point_id = 1
        point_map = Dict{Int, Int}()  # 原始点ID -> 新点ID
        
        # 写入原始晶粒顶点
        for (grain_id, vertex_ids) in grain_vertex_ids
            for id in vertex_ids
                if !haskey(point_map, id)
                    if haskey(original_points, id)
                        coord = original_points[id]
                        write(file, "Point($point_id) = {$(coord[1]), $(coord[2]), 0.0, lc};\n")
                        point_map[id] = point_id
                        point_id += 1
                    end
                end
            end
        end
        
        # 构建每个晶粒的顶点数组
        original_grain_points = Dict{Int, Vector{Int}}()
        for (grain_id, vertex_ids) in grain_vertex_ids
            original_grain_points[grain_id] = [point_map[id] for id in vertex_ids if haskey(point_map, id)]
        end
        
        # 写入偏移晶粒顶点
        write(file, "\n// 偏移点\n")
        offset_grain_points = Dict{Int, Vector{Int}}()
        
        for (grain_id, vertices) in offset_grains
            offset_points = Int[]
            for vertex in vertices
                write(file, "Point($point_id) = {$(vertex[1]), $(vertex[2]), 0.0, lc};\n")
                push!(offset_points, point_id)
                point_id += 1
            end
            offset_grain_points[grain_id] = offset_points
        end
        
        # ============ 写入线 ============
        write(file, "\n// 原始晶粒边界线\n")
        line_id = 1
        original_grain_lines = Dict{Int, Vector{Int}}()
        
        # 写入原始晶粒边界线
        for (grain_id, points) in original_grain_points
            n = length(points)
            if n < 3
                continue
            end
            
            grain_lines = Int[]
            for i in 1:n
                p1 = points[i]
                p2 = points[i % n + 1]
                write(file, "Line($line_id) = {$p1, $p2};\n")
                push!(grain_lines, line_id)
                line_id += 1
            end
            original_grain_lines[grain_id] = grain_lines
        end
        
        # 写入偏移晶粒边界线
        write(file, "\n// 偏移晶粒边界线\n")
        offset_grain_lines = Dict{Int, Vector{Int}}()
        
        for (grain_id, points) in offset_grain_points
            n = length(points)
            if n < 3
                continue
            end
            
            grain_lines = Int[]
            for i in 1:n
                p1 = points[i]
                p2 = points[i % n + 1]
                write(file, "Line($line_id) = {$p1, $p2};\n")
                push!(grain_lines, line_id)
                line_id += 1
            end
            offset_grain_lines[grain_id] = grain_lines
        end
        
        # 写入径向连接线（从老顶点到对应的新顶点）
        write(file, "\n// 老顶点到对应新顶点的径向连接线\n")
        radial_lines = Dict{Int, Dict{Int, Int}}()  # grain_id -> (orig_idx -> line_id)
        
        for (grain_id, orig_points) in original_grain_points
            if !haskey(offset_grain_points, grain_id) || !haskey(vertex_mappings, grain_id)
                continue
            end
            
            offset_points = offset_grain_points[grain_id]
            mapping = vertex_mappings[grain_id]
            
            grain_radial_lines = Dict{Int, Int}()
            
            for i in 1:length(orig_points)
                orig_pt = orig_points[i]
                
                # 使用映射找到对应的新顶点
                if haskey(mapping, i)
                    new_idx = mapping[i]
                    if 1 <= new_idx <= length(offset_points)
                        new_pt = offset_points[new_idx]
                        
                        write(file, "Line($line_id) = {$orig_pt, $new_pt}; // 顶点 $i 的径向连接\n")
                        grain_radial_lines[i] = line_id
                        line_id += 1
                    end
                end
            end
            
            radial_lines[grain_id] = grain_radial_lines
        end
        
        # ============ 写入线环 ============
        write(file, "\n// 线环\n")
        loop_id = 1
        
        # 偏移晶粒内部线环
        offset_loops = Dict{Int, Int}()
        
        for (grain_id, lines) in offset_grain_lines
            lines_str = join(lines, ", ")
            write(file, "Line Loop($loop_id) = {$lines_str};\n")
            offset_loops[grain_id] = loop_id
            loop_id += 1
        end
        
        # 晶界四边形线环
        boundary_loops = Dict{Tuple{Int, Int}, Int}()
        
        for (grain_id, orig_lines) in original_grain_lines
            if !haskey(offset_grain_lines, grain_id) || !haskey(radial_lines, grain_id)
                continue
            end
            
            offset_lines = offset_grain_lines[grain_id]
            rad_lines_map = radial_lines[grain_id]
            mapping = vertex_mappings[grain_id]
            
            n = length(orig_lines)
            
            for i in 1:n
                # 获取当前边的起点和终点索引
                start_idx = i
                end_idx = i % n + 1
                
                # 确保我们有这两个点的径向连接
                if haskey(rad_lines_map, start_idx) && haskey(rad_lines_map, end_idx)
                    orig_line = orig_lines[i]
                    
                    # 找到对应的偏移边
                    if haskey(mapping, start_idx) && haskey(mapping, end_idx)
                        offset_start_idx = mapping[start_idx]
                        offset_end_idx = mapping[end_idx]
                        
                        # 找到连接这两个偏移点的边
                        offset_line_idx = nothing
                        for j in 1:length(offset_lines)
                            if (j % length(offset_lines) + 1 == offset_end_idx && j == offset_start_idx) ||
                               (j == offset_end_idx && j % length(offset_lines) + 1 == offset_start_idx)
                                offset_line_idx = j
                                break
                            end
                        end
                        
                        if offset_line_idx !== nothing
                            offset_line = offset_lines[offset_line_idx]
                            rad_line1 = rad_lines_map[start_idx]
                            rad_line2 = rad_lines_map[end_idx]
                            
                            # 确保线的方向正确形成闭合四边形
                            # 顺序：原始边 -> 第二个顶点的径向连接 -> 偏移边(反向) -> 第一个顶点的径向连接(反向)
                            line_loop_str = "$orig_line, $rad_line2, -$offset_line, -$rad_line1"
                            write(file, "Line Loop($loop_id) = {$line_loop_str}; // 晶界四边形\n")
                            
                            boundary_loops[(grain_id, i)] = loop_id
                            loop_id += 1
                        end
                    end
                end
            end
        end
        
        # ============ 写入平面 ============
        write(file, "\n// 平面\n")
        surface_id = 1
        
        # 偏移晶粒内部平面
        offset_surfaces = Dict{Int, Int}()
        
        for (grain_id, loop_id) in offset_loops
            write(file, "Plane Surface($surface_id) = {$loop_id};\n")
            offset_surfaces[grain_id] = surface_id
            surface_id += 1
        end
        
        # 晶界四边形平面
        boundary_surfaces = Dict{Tuple{Int, Int}, Int}()
        
        for (key, loop_id) in boundary_loops
            write(file, "Plane Surface($surface_id) = {$loop_id};\n")
            boundary_surfaces[key] = surface_id
            surface_id += 1
        end
        
        # ============ 写入物理组 ============
        write(file, "\n// 物理组\n")
        
        # 所有偏移晶粒的物理组
        if !isempty(offset_surfaces)
            surfaces_str = join(values(offset_surfaces), ", ")
            write(file, "Physical Surface(\"AllGrains\") = {$surfaces_str};\n")
        end
        
        # 所有晶界四边形的物理组
        if !isempty(boundary_surfaces)
            surfaces_str = join(values(boundary_surfaces), ", ")
            write(file, "Physical Surface(\"AllBoundaries\") = {$surfaces_str};\n")
        end
        
        # 每个偏移晶粒的物理组
        for (grain_id, surface_id) in offset_surfaces
            write(file, "Physical Surface(\"Grain_$grain_id\") = {$surface_id};\n")
        end
        
        # 每个晶界四边形的物理组
        for (key, surface_id) in boundary_surfaces
            grain_id, edge_id = key
            write(file, "Physical Surface(\"Boundary_$(grain_id)_$edge_id\") = {$surface_id};\n")
        end
        
        # 径向连接线的物理组
        all_radial_lines = Int[]
        for grain_map in values(radial_lines)
            append!(all_radial_lines, values(grain_map))
        end
        
        if !isempty(all_radial_lines)
            lines_str = join(all_radial_lines, ", ")
            write(file, "Physical Line(\"RadialConnections\") = {$lines_str};\n")
        end
    end
    
    println("晶界模型已保存至: $filename")
end

"""
可视化晶界模型，突出显示老顶点与对应新顶点的连接关系
"""
function visualize_grain_boundaries(grains, offset_grains, vertex_mappings, output_file)
    p = plot(aspect_ratio=:equal, title="晶界模型 - 边平移算法", 
             xlabel="X", ylabel="Y", legend=true, size=(800, 800),
             dpi=300, grid=false, background_color=:white)
    
    # 绘制偏移晶粒
    for (grain_id, vertices) in offset_grains
        x = [v[1] for v in vertices]
        y = [v[2] for v in vertices]
        
        # 闭合多边形
        push!(x, x[1])
        push!(y, y[1])
        
        plot!(p, x, y, label=grain_id == first(keys(offset_grains)) ? "偏移晶粒" : "",
              color=:blue, linewidth=1.5)
    end
    
    # 绘制原始晶粒
    for (grain_id, vertices) in grains
        x = [v[1] for v in vertices]
        y = [v[2] for v in vertices]
        
        # 闭合多边形
        push!(x, x[1])
        push!(y, y[1])
        
        plot!(p, x, y, label=grain_id == first(keys(grains)) ? "原始晶粒" : "",
              color=:black, linewidth=1.5)
    end
    
    # 绘制径向连接线，突出显示老顶点与对应新顶点的连接
    for (grain_id, vertices) in grains
        if !haskey(offset_grains, grain_id) || !haskey(vertex_mappings, grain_id)
            continue
        end
        
        offset_vertices = offset_grains[grain_id]
        mapping = vertex_mappings[grain_id]
        
        for i in 1:length(vertices)
            if haskey(mapping, i)
                new_idx = mapping[i]
                if 1 <= new_idx <= length(offset_vertices)
                    orig = vertices[i]
                    offset = offset_vertices[new_idx]
                    
                    plot!(p, [orig[1], offset[1]], [orig[2], offset[2]],
                          label=i == 1 && grain_id == first(keys(grains)) ? "径向连接" : "",
                          color=:red, linewidth=1.5)
                end
            end
        end
    end
    
    savefig(p, output_file)
    println("可视化结果已保存至: $output_file")
end

"""
主函数
"""
function main()
    input_file = "equiaxed_grains.geo"
    output_file = "equiaxed_grains_with_boundaries.geo"
    
    # 检查输入文件是否存在
    if !isfile(input_file)
        println("错误: 输入文件 $input_file 不存在")
        println("当前目录: $(pwd())")
        println("目录内容: $(readdir())")
        return
    end
    
    # 偏移距离
    offset_distance = 0.002
    
    # 生成晶界模型
    generate_grain_boundary_model(input_file, output_file, offset_distance)
    
    println("处理完成！")
end

# 执行主函数
main()