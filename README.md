
# 多晶结构网格生成工具说明文档

本项目是一个用于生成多晶材料网格的工具集，基于Julia语言和Gmsh库实现。主要用于材料科学和计算力学领域的晶体结构模拟。
## 依赖包安装

本项目依赖以下 Julia 包：

1. **Gmsh** - 用于网格生成和操作
2. **Dates** - 用于日期和时间处理
3. **Statistics** - 用于统计计算
4. **DelimitedFiles** - 用于读写文本文件
5. **LinearAlgebra** - 用于矩阵运算
6. **Plots** (可选) - 用于可视化 (仅在某些脚本中使用)

### 安装方法

#### 方法一：使用项目的 Project.toml（推荐）

如果您已克隆此仓库，可以直接在项目目录中激活环境并安装所有依赖：

```julia
# 启动 Julia REPL，在项目目录中执行
julia

# 进入包管理模式（按 ]）
activate .
instantiate

```

#### 方法二：手动安装依赖包

或者，您也可以手动安装所需的包：

```julia
# 启动 Julia REPL
julia

# 进入包管理模式（按 ]）
add Gmsh
add Dates
add Statistics
add DelimitedFiles
add LinearAlgebra
# 可选：add Plots
# 按退格键返回 Julia REPL
```

### 关于 Gmsh

`Gmsh.jl` 包会自动下载并安装必要的 Gmsh 库文件，用户不需要单独安装 Gmsh 软件。但是，如果您想使用 Gmsh 的图形界面功能，或者想从命令行直接使用 Gmsh，您仍然可以从官方网站 https://gmsh.info/ 下载并安装完整的 Gmsh 软件。

## 文件概述

本项目包含以下几个主要的Julia脚本文件：

1. `mesh_grains_GB_alltri.jl` - 生成全三角形网格（晶界细密，晶粒粗糙）
2. `mesh_grains_GB_mix.jl` - 生成混合网格（晶界三角形，晶粒四边形和三角形混合）
3. `mesh_grains_GB_quadDOMI.jl` - 生成四边形主导网格（尽可能使用四边形，必要时使用三角形）
4. `Generate_seperateGBgroup.jl` - 生成分组的晶界，用于高级网格生成和控制

## 详细说明

### 1. mesh_grains_GB_alltri.jl

此文件借助输入的geo文件（由Generate_seperateGBgroup.jl生成）生成全三角形网格，特点是晶界区域网格细密，晶粒内部网格较粗糙。

**主要功能**：
- 生成纯三角形网格
- 晶界区域网格细化
- 晶粒内部网格粗化
- 两个区域之间平滑过渡

**使用方法**：
```julia
julia mesh_grains_GB_alltri.jl
```

### 2. mesh_grains_GB_mix.jl

此文件借助输入的geo文件（由Generate_seperateGBgroup.jl生成），生成混合类型网格，晶界区域使用三角形网格，晶粒内部使用四边形和三角形混合网格。

**主要功能**：
- 生成混合网格（三角形+四边形）
- 晶界区域使用三角形网格
- 晶粒内部尽可能使用四边形网格
- 复杂区域使用三角形网格填充

**使用方法**：
```julia
julia mesh_grains_GB_mix.jl
```

### 3. mesh_grains_GB_quadDOMI.jl

此文件借助输入的geo文件（由Generate_seperateGBgroup.jl生成），生成四边形主导网格，尽可能多地使用四边形网格单元，必要时使用三角形填充不规则区域。

**主要功能**：
- 生成四边形主导网格
- 使用高级拓扑算法识别适合四边形的区域
- 晶界和晶粒区域都尽可能使用四边形（晶界正常情况下全是四边形）
- 在不规则或复杂区域使用三角形填充

**使用方法**：
```julia
julia mesh_grains_GB_quadDOMI.jl
```

### 4. Generate_seperateGBgroup.jl

此文件用于产生有厚度晶界，并将信息记录进equiaxed_grains_with_boundaries.geo文件
需要输入geo文件，可以由neper等生成，或者自己编写。（需要在jl文件修改输入的文件名字）

**主要功能**：
- 根据边平移算法产生物理厚度的晶界
- 根据特定标准将晶界分组
- 生成包含分组信息的新几何geo文件

**使用方法**：
```julia
julia Generate_seperateGBgroup.jl
```

## 文件格式说明

本项目使用以下文件格式：

- `.geo` - Gmsh几何定义文件
- `.msh` - Gmsh网格文件
- `.inp` - Abaqus输入文件
- `.vtk` - VTK文件，用于后处理和可视化
- `.txt` - 文本信息文件，记录网格生成过程和统计信息

## 常见问题解答

1. **如何选择合适的网格类型？**
   - 三角形网格（alltri）：适合需要几何适应性高的情况，可以更好地拟合复杂几何形状
   - 混合网格（mix）：平衡了计算效率和几何适应性，一般情况下推荐使用
   - 四边形主导（quadDOMI）：适合需要高计算精度的情况，尤其是应变梯度较大的区域

2. **如何调整网格密度？**
   - 通过调整`grain_size`和`boundary_size`参数来控制不同区域的网格密度
   - `transition_distance`参数控制网格密度过渡区域的宽度

3. **如何处理大型模型？**
   - 增加`num_threads`参数来利用多核处理能力
   - 减少`max_optimization_time`以平衡网格质量和生成时间
   - 对于特别大的模型，可以考虑分区域生成然后合并

4. **如何导出不同格式的网格？**
   所有脚本都支持生成多种格式的网格文件，包括：
   - `.msh` - Gmsh原生格式
   - `.inp` - Abaqus输入文件
   - `.vtk` - 可视化工具包格式

## 使用建议

1. 首先使用`Generate_seperateGBgroup.jl`分析和准备几何模型
2. 根据具体需求选择合适的网格生成脚本
3. 有概率无法生成尺寸具有强烈对比的网格，请仔细选择参数。
4. 使用提供的`*_info.txt`文件检查网格质量和统计信息 


