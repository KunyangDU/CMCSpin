
using Clustering,MultivariateStats
include("../src/src.jl")

function FCM(data::Matrix, centers::Matrix, m=2.0)
    n_samples = size(data, 2)
    n_clusters = size(centers, 2)
    weights = zeros(n_clusters, n_samples)
    
    for i in 1:n_samples
        x = data[:, i]
        dists = [norm(x - centers[:, c]) for c in 1:n_clusters]
        
        # FCM 隶属度公式
        for c in 1:n_clusters
            # 避免除以 0
            if dists[c] < 1e-10
                weights[:, i] .= 0.0
                weights[c, i] = 1.0
                break
            end
            
            sum_denom = sum((dists[c] / dists[j])^(2/(m-1)) for j in 1:n_clusters)
            weights[c, i] = 1.0 / sum_denom
        end
    end
    return weights
end

"""
    find_peaks_simple(y::AbstractVector; min_height=0.0, min_dist=1)

在不引入新包的情况下寻找峰值。
策略：找到所有局部极大值 -> 按高度降序排列 -> 剔除距离过近的次要峰。

参数:
- `y`: 输入的一维数据（假设已平滑）
- `min_height`: 峰值的最小绝对高度（用于过滤背景/噪声）
- `min_dist`: 两个峰之间的最小索引距离（防止对同一个宽峰重复计数）

返回:
- `peaks`: 峰值在 y 中的索引列表 (按索引从小到大排序)
"""
function find_peaks_simple(y::AbstractVector{T}; min_height=0.0, min_dist=1) where T<:Real
    len = length(y)
    if len < 3
        return Int[]
    end

    # 1. 第一轮扫描：寻找简单的局部极大值候选点
    candidates = Int[]
    
    # 简单的三点比较 (i-1 < i >= i+1)
    # 使用 >= i+1 是为了处理平顶 (Plateau) 的情况，取平顶的左边缘
    for i in 2:len-1
        if y[i] > y[i-1] && y[i] >= y[i+1]
            if y[i] > min_height
                push!(candidates, i)
            end
        end
    end

    # 如果没有找到或者不需要距离过滤，直接返回
    if isempty(candidates) || min_dist <= 1
        return candidates
    end

    # 2. 第二轮筛选：基于距离的贪婪剔除 (Greedy Suppression)
    
    # 将候选点按高度从大到小排序
    # 这样保证我们优先保留最高的峰，剔除它旁边矮的假峰
    sort!(candidates, by = i -> y[i], rev=true)

    kept_peaks = Int[]
    
    for cand in candidates
        is_isolated = true
        for existing in kept_peaks
            # 如果当前候选点和已经保留的高峰距离太近，就丢弃它
            if abs(cand - existing) < min_dist
                is_isolated = false
                break
            end
        end
        
        if is_isolated
            push!(kept_peaks, cand)
        end
    end

    # 3. 最后按索引顺序（即温度顺序）返回
    return sort!(kept_peaks)
end


T0 = 2.0
Tf = 0.5
α = 0.9
params = (
    J = 1.0,
)

θ = 0.0 * pi
ϕ = 0.0 * pi

lsHf = 0.0:0.2:8.0
SSF = zeros(length(Latt),length(lsHf))
for (i,Hf) in enumerate(lsHf)
# Hf = 8.0
@load "lab/data/ψ_$(round(θ/pi))_$(round(ϕ/pi))_$(Hf)_$(params)_$(T0)_$(Tf)_$(α).jld2" ψ

SM = ψ.pattern' * ψ.pattern
lspq = [[nx//Lx,ny//Ly] for nx in 0:Lx-1 for ny in 0:Ly-1]
lsk = 2pi .* lspq
lstk = Tuple.(lsk)
SSF[:,i] = FT2(Latt,SM,lstk)
end

pca_model = fit(PCA, SSF;maxoutdim = 5)
SSF_proj = MultivariateStats.transform(pca_model, SSF)
fd = [dot(SSF[:,i],SSF[:,i+1]) for i in 1:size(SSF)[2]-1]
W = (SSF' * SSF) .^ 4
D = diagm(sum(W,dims = 1)[:])
f = eigen(D-W)
K = findmax(diff(f.values)[1:10])[2]

# maximum(diff(f.values[2:end]))
# minimum(W)
# fv = norm.(eachcol(diff(SSF,dims = 2)))[:]
# fd = [(1 - dot(SSF[:,i],SSF[:,i+1])^2)/(lsHf[i+1] - lsHf[i])^2 for i in 1:size(SSF)[2]-1]
# pinds = find_peaks_simple(fd)


# clustering = kmeans(SSF_proj,k)
# centers = clustering.centers
# weights = FCM(SSF_proj,clustering.centers)

# assignments(clustering)
# FCM_projection(SSF,clustering.centers)
# clustering = fuzzy_cmeans(SSF,2,2.0)
# clustering.centers
# clustering.weights

# clustering = dbscan(SSF, 0.07, min_neighbors = 2)

# clustering.clusters
# SSF
# pca_model = fit(PCA, reshape(SSF,1,:))
# SSF_proj = MultivariateStats.transform(pca_model, reshape(SSF,1,:))
# SSF
# fig = Figure()
# ax = Axis(fig[1,1])
# x = map(x -> x[1],lsk)
# y = map(x -> x[2],lsk)

# heatmap!(ax, x,y,SSF)
# display(fig)

fig = Figure()
ax = Axis(fig[1,1])
scatterlines!(ax,f.values[1:10])
display(fig)
# fv
# dot(SSF[:,1],SSF[:,2])
# fv
# weights