using FFTW
using LinearAlgebra

# 1. 定义晶格参数
Nx = 6  # x方向格点数
Ny = 6  # y方向格点数

# 2. 定义实空间函数 f(r)
# 这里作为一个示例，我们定义一个高斯波包或者简单的余弦波
# r = (x, y) 索引从 1 到 N
real_space_data = zeros(ComplexF64, Nx, Ny)
for x in 1:Nx, y in 1:Ny
    # 示例：一个简单的波包
    real_space_data[x, y] = exp(-((x - Nx/2)^2 + (y - Ny/2)^2) / 4.0)
end

# 3. 进行快速傅里叶变换 (FFT)
# 注意：Julia 的 fft 默认没有归一化因子 (如 1/sqrt(N))，逆变换 ifft 有 1/N 因子。
# 物理中通常为了保范数(Unitary)，手动乘以 1/sqrt(Nx*Ny)
k_space_data = fft(real_space_data) ./ sqrt(Nx * Ny)

# ==========================================
# 关键部分：构造有效的 q 点 (Reciprocal Lattice Vectors)
# ==========================================

# 方法 A: 标准 FFT 顺序 (Standard Order)
# 对应 fft 输出的索引：[0, 1, ..., N/2-1, -N/2, ..., -1] * (2π/N)
# 范围大致是 [0, 2π)
qx_std = 2π * fftfreq(Nx)
qy_std = 2π * fftfreq(Ny)
k_space_data
# println("--- 标准 FFT 顺序 qx (前几个) ---")
# println(qx_std[1:3]) 

# # 方法 B: 中心化顺序 (Shifted Order) - 物理常用
# # 将零频 (q=0) 移到数组中心。范围大致是 [-π, π)
# # 这对于绘制能带图或结构因子 S(q) 非常有用。
# qx_shifted = 2π * fftshift(fftfreq(Nx))
# qy_shifted = 2π * fftshift(fftfreq(Ny))
# k_space_data_shifted = fftshift(k_space_data)

# println("\n--- 中心化 (Shifted) qx ---")
# println(qx_shifted)

# # 4. 如何访问特定 q 点的数据
# # 生成 q 点的网格 (Meshgrid)
# Q_grid = [(qx, qy) for qx in qx_shifted, qy in qy_shifted]

# # 示例：打印 q=(0,0) 处的傅里叶分量幅值
# # 在 shifted 数组中，中心点索引为:
# cx = div(Nx, 2) + 1
# cy = div(Ny, 2) + 1
# println("\n--- q=(0,0) 处的分量 ---")
# println("q vector: ", Q_grid[cx, cy])
# println("Amplitude: ", abs(k_space_data_shifted[cx, cy]))