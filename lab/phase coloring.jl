using CairoMakie
using Colors
using ColorSchemes

# 1. 模拟数据：假设有 5 个相 (L x L x 5)
L = 100
N_phases = 5
weights = rand(L, L, N_phases)

# 物理约束模拟：
# 通常在某一点主要是一个相占主导，或者两个相竞争
# 这里为了让图好看，做一次 softmax 或者归一化，模拟 Order Parameter
weights .= weights .^ 3 # 锐化一下，模拟相分界
for i in 1:L, j in 1:L
    weights[i,j,:] ./= sum(weights[i,j,:])
end

# ---------------------------------------------------------
# 2. 定义基矢颜色 (Basis Colors)
# ---------------------------------------------------------
# 使用 ColorSchemes 中区分度高的离散色盘，例如 tab10, Dark2, Paired
# 也可以手动指定 colors = [colorant"red", colorant"blue", ...]
basis_colors = ColorSchemes.tab10[1:N_phases]

# ---------------------------------------------------------
# 3. 颜色混合核心逻辑 (Linear Combination)
# ---------------------------------------------------------
# 将权重张量 (L,L,N) 压缩为颜色矩阵 (L,L)
# 计算公式：RGB_final = sum( weight_k * RGB_basis_k )

function mix_phases(weights, bases)
    (L_x, L_y, N) = size(weights)
    img = Matrix{RGB{Float64}}(undef, L_x, L_y)
    
    # 预先提取基矢的分量，加速循环
    rs = red.(bases)
    gs = green.(bases)
    bs = blue.(bases)
    
    for i in 1:L_x, j in 1:L_y
        # 标量积：w · r_basis
        # 这里的 view(weights, i, j, :) 避免了内存分配，或者直接写循环
        w = view(weights, i, j, :)
        
        r = sum(w .* rs)
        g = sum(w .* gs)
        b = sum(w .* bs)
        
        # 截断以防数值误差导致溢出 [0,1]
        img[i,j] = RGB(clamp(r,0,1), clamp(g,0,1), clamp(b,0,1))
    end
    return img
end

img_matrix = mix_phases(weights, basis_colors)

# ---------------------------------------------------------
# 4. 绘图与图例
# ---------------------------------------------------------
fig = Figure(size = (800, 600))

# 4.1 绘制相图
ax = Axis(fig[1, 1], 
    title = "Multi-Phase Diagram (N=$N_phases)",
    aspect = DataAspect(),
    xlabel = "Parameter 1", ylabel = "Parameter 2"
)

heatmap!(ax, 1:L, 1:L, img_matrix, interpolate = false)

# 4.2 构造图例 (Legend)
# 对于多相系统，图例至关重要
legend_elements = [PolyElement(color = c, polygon = Point2f[(0,0), (1,0), (1,1), (0,1)]) for c in basis_colors]
legend_labels = ["Phase $k" for k in 1:N_phases]

Legend(fig[1, 2], legend_elements, legend_labels, "Phase Bases")

display(fig)