using LinearAlgebra
using CairoMakie

# ==============================================================================
# 1. 核心求解器：给定构型，计算色散
# ==============================================================================

"""
    Configuration
结构体用于存储晶格信息和自旋构型
"""
struct SpinConfig
    Lx::Int
    Ly::Int
    S_val::Float64
    J::Float64
    spins::Matrix{Vector{Float64}} # 维度 (Lx, Ly)，每个元素是 [sx, sy, sz]
end

"""
    get_local_basis(n)
根据自旋方向 n (z_local)，构建局部坐标系 (x_local, y_local, z_local)
"""
function get_local_basis(n::Vector{Float64})
    z_loc = n / norm(n)
    # 选取辅助向量构建 x_loc，防止与 z_loc 平行
    aux = abs(z_loc[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    x_loc = cross(aux, z_loc)
    x_loc /= norm(x_loc)
    y_loc = cross(z_loc, x_loc)
    return (x_loc, y_loc, z_loc)
end

"""
    calculate_lswt_dispersion(config, k_path)
输入:
- config: 自旋构型结构体
- k_path: 一系列 k 点坐标 [[kx, ky], ...]
输出:
- freqs: 对应每个 k 点的所有能级 (Sorted)
"""
function calculate_lswt_dispersion(cfg::SpinConfig, k_path::Vector{Vector{Float64}})
    Lx, Ly = cfg.Lx, cfg.Ly
    N_sites = Lx * Ly
    S = cfg.S_val
    J = cfg.J
    
    # 1. 预计算每个格点的局部坐标系
    # basis[i, j] = (x_hat, y_hat, z_hat)
    basis = [get_local_basis(cfg.spins[i, j]) for i in 1:Lx, j in 1:Ly]

    # 将 2D 坐标 (ix, iy) 映射到 1D 索引 1..N
    # Julia 是列优先 (Column-major)，所以 index = (iy-1)*Lx + ix
    to_idx(ix, iy) = (iy - 1) * Lx + ix

    n_k = length(k_path)
    all_energies = zeros(Float64, n_k, N_sites) # 只取正能级，共 N 个分支

    # 定义最近邻偏移 (正方晶格)
    deltas = [(1, 0), (0, 1), (-1, 0), (0, -1)]

    for (k_idx, k_vec) in enumerate(k_path)
        kx, ky = k_vec
        
        # 构造 Bogoliubov 矩阵 M = 2N x 2N
        # 块结构: [ A   B ]
        #         [ B* A*]
        A = zeros(ComplexF64, N_sites, N_sites)
        B = zeros(ComplexF64, N_sites, N_sites)
        
        # 遍历原胞内所有格点
        for ix in 1:Lx, iy in 1:Ly
            u = to_idx(ix, iy)
            x_u, y_u, z_u = basis[ix, iy] # Site u 的局部基矢
            
            # 遍历邻居
            for (dx, dy) in deltas
                # 目标格点坐标（处理周期性边界）
                jx_raw = ix + dx
                jy_raw = iy + dy
                
                # 确定相位因子 (Bloch Phase)
                # 如果跳出了原胞，不仅要取模找到对应的 j，还要乘上 exp(i k R)
                phase = 0.0
                
                # 处理 X 方向边界
                jx = mod1(jx_raw, Lx)
                if jx_raw > Lx; phase += kx * Lx; end
                if jx_raw < 1;  phase -= kx * Lx; end
                
                # 处理 Y 方向边界
                jy = mod1(jy_raw, Ly)
                if jy_raw > Ly; phase += ky * Ly; end
                if jy_raw < 1;  phase -= ky * Ly; end
                
                v = to_idx(jx, jy)
                x_v, y_v, z_v = basis[jx, jy] # Site v 的局部基矢
                
                # 相因子 e^{-i k (R_u - R_v)}? 
                # 标准写法：算符傅里叶变换后，hopping term 携带 e^{i k (r_j - r_i)}
                # 这里 r_j - r_i 正好是跨越原胞的位移向量
                phase_factor = exp(im * phase) 
                
                # === 计算相互作用矩阵元 ===
                # Hamiltonian term: J * S_u . S_v
                # 我们需要计算 S_u^+ S_v^-, S_u^+ S_v^+ 等的系数
                
                # 为了通用，直接利用基矢点积展开：
                # S_u . S_v = (x_u S_u^x + y_u S_u^y + z_u S_u^z) . (x_v ... )
                # 其中 S^z -> S - a'a, S^x -> sqrt(S/2)(a+a'), S^y -> sqrt(S/2)/i (a-a')
                
                # 1. 对角项 contributions (from S_u^z S_v^z)
                # S_u^z S_v^z ~= S^2 - S(n_u + n_v)
                # 贡献到 A[u,u] 和 A[v,v] (注意求和重复计算的问题，这里我们遍历每条键的一端)
                # 这里的循环是 "对每个u，遍历其邻居v"，所以是定向键。系数取 J/2 * 2 = J
                
                zz_dot = dot(z_u, z_v)
                diag_term = -0.5 * J * S * zz_dot 
                # Explanation: H_linearized has term -J*S * (z_u.z_v) * (a_u^\dagger a_u)
                A[u, u] += diag_term
                
                # 2. Hopping (A matrix): a_u^\dagger a_v
                # 来自 S_u^x S_v^x, S_u^y S_v^y 等混合项
                # Coeff = (J*S/2) * [ (x_u - i y_u) . (x_v + i y_v) ] * phase
                # 对应算符 a_u^\dagger a_v
                
                vec_u_minus = x_u - im * y_u # for a_u^\dagger
                vec_v_plus  = x_v + im * y_v # for a_v
                
                coeff_hop = (0.25 * J * S) * sum(vec_u_minus .* vec_v_plus) * phase_factor
                A[u, v] += coeff_hop
                
                # 3. Pairing (B matrix): a_u a_v
                # Coeff = (J*S/2) * [ (x_u + i y_u) . (x_v + i y_v) ] * phase
                # 对应算符 a_u a_v
                # 注意：最终矩阵形式是 1/2 [a' a] H [a; a']，所以这里写入 B[u,v] 
                
                vec_u_plus = x_u + im * y_u
                # vec_v_plus defined above
                
                coeff_pair = (0.25 * J * S) * sum(vec_u_plus .* vec_v_plus) * phase_factor
                B[u, v] += coeff_pair
            end
        end
        
        # 构造 BdG 矩阵并对角化
        # M = g * H_total
        # H_total = [A B; B' A'] (Hermitian)
        # g = diag(1, -1)
        H_total = [A B; adjoint(B) conj(A)]
        g = Diagonal([ones(N_sites); -ones(N_sites)])
        
        M = g * H_total
        
        # 求解非厄米本征值
        evals = eigen(M).values
        
        # 筛选实部为正的本征值 (物理激发)
        # 稍微加一点容差处理数值误差
        real_energies = sort(real.(filter(e -> real(e) > 1e-6, evals)))
        
        # 存储前 N_sites 个模式（通常就是全部正能级）
        if length(real_energies) >= N_sites
            all_energies[k_idx, :] = real_energies[1:N_sites]
        end
    end
    
    return all_energies
end

# ==============================================================================
# 2. 示例：Néel 反铁磁 (AFM) 构型验证
# ==============================================================================

# 1. 设定参数
Lx, Ly = 1, 1  # 磁性原胞大小 (AFM 至少需要 2x2 或 2x1)
J = 1.0
S = 0.5

# 2. 构造 AFM 自旋构型 (Up-Down Checkerboard)
spins = Matrix{Vector{Float64}}(undef, Lx, Ly)
for i in 1:Lx, j in 1:Ly
    # (-1)^(i+j) 决定方向
    if iseven(i + j)
        spins[i, j] = [0.0, 0.0, 1.0]  # Up
    else
        spins[i, j] = [0.0, 0.0, 1.0] # Down
    end
end

config = SpinConfig(Lx, Ly, S, J, spins)

# 3. 定义 K 点路径 (High symmetry points for Square Lattice)
# Gamma -> X -> M -> Gamma
nk = 100
path1 = [[k, 0.0] for k in range(0, π, length=nk)]        # G -> X
path2 = [[π, k] for k in range(0, π, length=nk)]          # X -> M
path3 = [[k, k] for k in range(π, 0, length=nk)]          # M -> G

k_path_vec = vcat(path1, path2, path3)
x_axis = 1:length(k_path_vec)

# 4. 计算
bands = calculate_lswt_dispersion(config, k_path_vec)

fig = Figure()
ax = Axis(fig[1,1])
for bd in eachcol(bands)
    lines!(ax, bd)
end
display(fig)

# # 5. 绘图
# plot(x_axis, bands, 
#      label=["Band 1" "Band 2" "Band 3" "Band 4"], 
#      linewidth=2, 
#      title="Spin Wave Dispersion (Square Lattice AFM)",
#      xlabel="Momentum Path (G-X-M-G)", ylabel="Energy (J)",
#      legend=:topright)

# # 标记高对称点
# vline!([nk, 2nk], label="", color=:black, linestyle=:dash)
# xticks!([1, nk, 2nk, 3nk], ["Γ", "X", "M", "Γ"])