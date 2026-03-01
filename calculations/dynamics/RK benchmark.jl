using LinearAlgebra
# using Printf

# 1. 定义物理方程: dS/dt = S x B
# 为了简单，直接把 B 写死在函数里 [0, 0, 1]
function dynamics(S)
    B = [0.0, 0.0, 1.0]
    return cross(S, B) # 叉乘
end

# 2. 标准 RK4 单步函数
function rk4(y::Union{Number, Vector}, dt::Number, f::Function)
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    
    return @. y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
end

# 3. 主程序
function main()
    # 初始条件
    S = [1.0, 0.0, 0.0]  # 初始自旋在 x 轴
    dt = 0.01            # 步长
    t_end = 10.0         # 总时间
    steps = Int(t_end / dt)

    println("开始积分...")
    println("初始状态: $S")
    
    # 时间循环
    for i in 1:steps
        S = step_rk4(S, dt, dynamics)
    end

    println("积分完成 (t = $t_end)")
    println("最终结果: $(S)\n")
    
    # 理论值对比 (绕 Z 轴旋转了 10 弧度，方向为 -10)
    # S_exact = [cos(-10), sin(-10), 0]
    s_exact_x = cos(-10.0)
    s_exact_y = sin(-10.0)
    println("---")
    println("理论真值: $([s_exact_x, s_exact_y, 0.0])\n")
end

main()