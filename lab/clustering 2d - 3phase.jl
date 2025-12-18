using LinearAlgebra,MultivariateStats,Clustering,Colors
include("../src/src.jl")

# function mix_phases(weights::Matrix, bases::Vector{RGB{Float64}} = ColorSchemes.tab10[1:size(weights)[1]])
#     (N,L) = size(weights)
#     img = Vector{RGB{Float64}}(undef, L)
    
#     # 预先提取基矢的分量，加速循环
#     rs = red.(bases)
#     gs = green.(bases)
#     bs = blue.(bases)
    
#     for i in 1:L
#         # 标量积：w · r_basis
#         # 这里的 view(weights, i, j, :) 避免了内存分配，或者直接写循环
#         w = weights[:,i]
        
#         r = sum(w .* rs)
#         g = sum(w .* gs)
#         b = sum(w .* bs)
        
#         # 截断以防数值误差导致溢出 [0,1]
#         img[i] = RGB(clamp(r,0,1), clamp(g,0,1), clamp(b,0,1))
#     end
#     return img
# end

function mix_phases(weights::Array{T,3}, bases::Vector = ColorSchemes.tab10[1:size(weights)[1]]) where T
    (N,L_x, L_y) = size(weights)
    img = Matrix{RGB{Float64}}(undef, L_x, L_y)
    
    # 预先提取基矢的分量，加速循环
    rs = red.(bases)
    gs = green.(bases)
    bs = blue.(bases)
    
    for i in 1:L_x, j in 1:L_y
        # 标量积：w · r_basis
        # 这里的 view(weights, i, j, :) 避免了内存分配，或者直接写循环
        w = view(weights, :, i, j)
        
        r = sum(w .* rs)
        g = sum(w .* gs)
        b = sum(w .* bs)
        
        # 截断以防数值误差导致溢出 [0,1]
        img[i,j] = RGB(clamp(r,0,1), clamp(g,0,1), clamp(b,0,1))
    end
    return img
end
# ===========================
# 1. 参数设置
# ===========================
L = 60           # 扫描尺寸 20x20
D = L^2        # 矢量维度 (模拟结构因子的大小，比如 q 点的数量)
width = 1.0      # Crossover 的宽度 (类似温度)。
                 # 调大(如 5.0)则混合区变大；调小(如 0.5)则边界锐利。

# ===========================
# 2. 定义三个纯相的基底向量 (相互正交)
# ===========================
# 模拟三个完全不同的序，比如 (pi,pi), (0,0), (pi,0)
# 这里简单起见，让它们在不同维度上有峰值
basis_1 = zeros(D); basis_1[10:20] .= 1.0; (basis_1) |> x -> x/norm(x)
basis_2 = zeros(D); basis_2[40:50] .= 1.0; (basis_2) |> x -> x/norm(x)
basis_3 = zeros(D); basis_3[70:80] .= 1.0; (basis_3) |> x -> x/norm(x)

bases = [basis_1, basis_2, basis_3]

# ===========================
# 3. 定义三个相在相图上的中心位置
# ===========================
# 设想呈三角形分布
centers = [
    (L/4, 0.0),    # 相 1 在左下
    (L, L/2),   # 相 2 在右下
    (0.0, L)   # 相 3 在上方
]

# ===========================
# 4. 生成数据
# ===========================
# data_grid[i, j] 存储一个长度为 D 的向量
data_grid = Matrix{Vector{Float64}}(undef, L, L)
# label_map 用于画图验证 (存 RGB 颜色)
rgb_map = zeros(L, L, 3) 

for i in 1:L
    for j in 1:L
        # 计算当前点 (i,j) 到三个中心的欧几里得距离
        dists = [norm([i, j] .- c) for c in centers]
        
        # 核心逻辑：使用 Softmax 机制产生 Crossover
        # 距离越近，权重(weight)越大
        weights = exp.(-dists ./ width)
        weights /= sum(weights) # 归一化权重
        
        # 混合向量：线性叠加 (模拟量子态叠加或经典统计平均)
        vec_ij = weights[1] * basis_1 + weights[2] * basis_2 + weights[3] * basis_3
        
        # 归一化 (通常结构因子或自旋都有归一化约束)
        (vec_ij) |> x -> x/norm(x)
        
        data_grid[i, j] = vec_ij
        
        # 记录权重用于可视化 (R, G, B 对应三个相的成分)
        rgb_map[j, i, :] = weights # 注意 heatmap 坐标系转置
    end
end


data = hcat(data_grid[:]...)


# pca_model = fit(PCA, data)
# data_proj = MultivariateStats.transform(pca_model, data)

# @time W = (data' * data) .^ 3
# @time D = diagm(sum(W,dims = 1)[:])
# # @time vals, vecs, info = eigsolve(D-W,ones(length(data_grid)),5,:SR)
# @time f = eigvals(D-W)
# K = findmax(diff(f)[1:5])[2]

pca_model = fit(PCA, data)
data_proj = MultivariateStats.transform(pca_model, data)
K = length(pca_model.prinvars) + 1
clustering = kmeans(data,K)
centers = clustering.centers
weights = FCM(data,clustering.centers)
weights = reshape(weights,3,L,L)
colors = mix_phases(weights,[colorant"red",colorant"green",colorant"blue"])
# weights[:,1,1]
# d = 3
coor = [[i,j] for i in 1:L for j in 1:L]
x = map(x -> x[1], coor)
y = map(x -> x[2], coor)
# colors = mix_phases(weights)
fig = Figure()
# heatmap!(Axis(fig[1,1]),1:L,1:L,rgb_map[:,:,d])
# heatmap!(Axis(fig[1,1]),1:L,1:L,colors')
scatter!(Axis(fig[1,1]), x,y,color = [:red,:blue,:green][clustering.assignments])
# scatterlines!(ax,f[1:5])
display(fig)
# length(x)

# clustering.

