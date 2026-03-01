using LinearAlgebra
using CairoMakie
using Statistics

# ==========================================
# 1. 物理模型部分 (与之前保持一致)
# ==========================================
struct SpinConfig
    Lx::Int
    Ly::Int
    S_val::Float64
    J::Float64
    spins::Matrix{Vector{Float64}}
end

# 辅助：局部坐标基矢
function get_local_basis(n::Vector{Float64})
    z_loc = n / norm(n)
    aux = abs(z_loc[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    x_loc = cross(aux, z_loc); x_loc /= norm(x_loc)
    y_loc = cross(z_loc, x_loc)
    return (x_loc, y_loc, z_loc)
end

"""
计算原始的离散能谱数据。
返回: (k_indices, energies, intensities) 的列表，用于后续展宽。
"""
function calculate_raw_spectrum(cfg::SpinConfig, k_path::Vector{Vector{Float64}})
    Lx, Ly = cfg.Lx, cfg.Ly
    N_sites = Lx * Ly
    S, J = cfg.S_val, cfg.J
    
    basis = [get_local_basis(cfg.spins[i, j]) for i in 1:Lx, j in 1:Ly]
    to_idx(ix, iy) = (iy - 1) * Lx + ix
    
    # 构建实空间 Ham (超胞 Gamma 点)
    A = zeros(ComplexF64, N_sites, N_sites)
    B = zeros(ComplexF64, N_sites, N_sites)
    deltas = [(1, 0), (0, 1), (-1, 0), (0, -1)]

    for ix in 1:Lx, iy in 1:Ly
        u = to_idx(ix, iy)
        x_u, y_u, z_u = basis[ix, iy]
        
for (dx, dy) in deltas
            jx = mod1(ix + dx, Lx); jy = mod1(iy + dy, Ly)
            v = to_idx(jx, jy)
            x_v, y_v, z_v = basis[jx, jy]
            
            # --- 修正开始 ---
            
            # 1. 对角项 (S^z S^z)
            zz = dot(z_u, z_v) # 实向量 dot 没问题
            A[u, u] += -0.5 * J * S * zz
            
            # 2. 构造复矢量 (不取共轭)
            # e_minus = x - iy (对应 S-)
            # e_plus  = x + iy (对应 S+)
            vec_u_minus = x_u - im * y_u
            vec_u_plus  = x_u + im * y_u
            
            vec_v_minus = x_v - im * y_v
            vec_v_plus  = x_v + im * y_v
            
            # 3. 计算 A 矩阵 (Hopping: a_u^dag a_v)
            # 对应算符项: S_u^- S_v^+ 
            # 几何系数: e_u^- \cdot e_v^+
            # 必须使用 sum( .* ) 避免 Julia 自动共轭第一个参数
            term_A = sum(vec_u_minus .* vec_v_plus)
            A[u, v] += (0.25 * J * S) * term_A
            
            # 4. 计算 B 矩阵 (Pairing: a_u a_v)
            # 对应算符项: S_u^+ S_v^+
            # 几何系数: e_u^+ \cdot e_v^+
            term_B = sum(vec_u_plus .* vec_v_plus)
            B[u, v] += (0.25 * J * S) * term_B
            
            # --- 修正结束 ---
        end
    end

    H_total = [A B; adjoint(B) conj(A)]
    g = Diagonal([ones(N_sites); -ones(N_sites)])
    M = g * H_total
    evals, V = eigen(M)
    
    # 筛选物理模式
    phys_mask = real(evals) .> 1e-4
    energies_raw = real(evals[phys_mask])
    vectors = V[:, phys_mask]
    
    # 结果容器
    k_list = Int[]
    E_list = Float64[]
    I_list = Float64[]
    
    r_vecs = [[ix, iy] for iy in 1:Ly for ix in 1:Lx]

    # 遍历 K 路径计算强度
    for (k_idx, Q_vec) in enumerate(k_path)
        Qx, Qy = Q_vec
        
        for n in 1:length(energies_raw)
            E = energies_raw[n]
            u_vec = vectors[1:N_sites, n]
            v_vec = vectors[N_sites+1:end, n]
            
            # 计算结构因子
            ft_sum = ComplexF64(0.0)
            for i in 1:N_sites
                ix, iy = r_vecs[i]
                phase = exp(-im * (Qx * ix + Qy * iy))
                # 简化近似: 假设测量垂直分量
                ft_sum += phase * (u_vec[i] + v_vec[i])
            end
            intensity = abs2(ft_sum)
            
            # 只有强度显著的点才记录，节省计算
            if intensity > 1e-3
                push!(k_list, k_idx)
                push!(E_list, E)
                push!(I_list, intensity)
            end
        end
    end
    
    return k_list, E_list, I_list
end

# ==========================================
# 2. 展宽处理：生成均匀网格数据 S(k, w)
# ==========================================
"""
将离散的 (k, E, I) 数据展宽到均匀的 (w_axis) 网格上
"""
function broaden_spectrum(k_indices, E_data, I_data, nk, w_axis, eta=0.1)
    nw = length(w_axis)
    # 创建矩阵 (nk, nw)，这是 Makie heatmap 喜欢的格式
    # S[k, w]
    S_matrix = zeros(Float64, nk, nw)
    
    # 洛伦兹函数 L(x) = (1/pi) * (eta / (x^2 + eta^2))
    # 为了速度，只在能量附近的窗口内计算
    
    for (k, E0, I0) in zip(k_indices, E_data, I_data)
        # 找到能量在 w_axis 上的大致位置
        # 简单暴力法：遍历所有 w (如果 w 轴很长，可以优化范围)
        for (w_idx, w) in enumerate(w_axis)
            diff = w - E0
            # 截断：太远的地方设为0，加速计算
            if abs(diff) < 10 * eta 
                val = (1/π) * (eta / (diff^2 + eta^2))
                S_matrix[k, w_idx] += I0 * val
            end
        end
    end
    
    return S_matrix
end

# ==========================================
# 3. 主程序与 CairoMakie 绘图
# ==========================================

# --- 准备数据 ---
Lx, Ly = 50, 10
spins = Matrix{Vector{Float64}}(undef, Lx, Ly)
for i in 1:Lx, j in 1:Ly
    spins[i, j] = iseven(i+j) ? [0.0,0.0,1.0] : [0.0,0.0,-1.0]
end
# 加上一点随机扰动模拟 MC 数据的不完美性
# for i in 1:Lx, j in 1:Ly; spins[i,j] += 0.05*randn(3); spins[i,j]/=norm(spins[i,j]); end

config = SpinConfig(Lx, Ly, 0.5, -1.0, spins)

nk_seg = 100
path_vec = vcat(
    # [[k, 0.0] for k in range(0, π, length=nk_seg)],
    [[k, pi] for k in range(0, π, length=nk_seg)],
    # [[π, k] for k in range(0, π, length=nk_seg)],
    # [[k, k] for k in range(π, 0, length=nk_seg)]
)
total_nk = length(path_vec)

println("1. Calculating eigenvalues...")
k_idxs, Es, Is = calculate_raw_spectrum(config, path_vec)

println("2. Broadening spectrum...")
# 定义能量轴 (Uniform Grid)
w_max = 2.5 # 根据 J 和 S 估算，AFM 最大能量通常是 ~4JS
w_axis = range(0, w_max, length=400)
eta = 0.08 # 展宽宽度

S_kw = broaden_spectrum(k_idxs, Es, Is, total_nk, w_axis, eta) / Lx / Ly

println("3. Plotting with CairoMakie...")

# --- CairoMakie 绘图 ---
fig = Figure(size = (600, 400), fontsize = 18)

# 创建坐标轴
ax = Axis(fig[1, 1],
    title = "Dynamical Structure Factor S(k, ω)",
    xlabel = "Momentum Path",
    ylabel = "Energy (J)",
    xticks = ([1, nk_seg, 2*nk_seg, 3*nk_seg], ["Γ", "X", "M", "Γ"])
)

# 绘制 Heatmap
# 注意：heatmap! 默认 x 对应第一维，y 对应第二维
# S_kw 是 (nk, nw)，这里 x=1:nk, y=w_axis
hm = heatmap!(ax, 1:total_nk, w_axis, S_kw,colorrange = (0,4))

# 添加垂直虚线标记高对称点
vlines!(ax, [nk_seg, 2*nk_seg], color = :white, linestyle = :dash, linewidth=1.5)

# 添加 Colorbar
Colorbar(fig[1, 2], hm, label = "Intensity (a.u.)")

# 显示或保存
display(fig) 
save("lswt_dispersion_heatmap_afm.png", fig)
# println("Done! Saved to lswt_dispersion_heatmap.png")