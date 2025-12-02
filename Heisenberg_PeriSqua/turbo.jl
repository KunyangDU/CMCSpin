using LoopVectorization

# 准备数据：1000万个点
N = 1000000000
a = rand(N)
b = rand(N)
y = zeros(N)

# 1. 原生 Julia 写法
function heavy_calc_base!(y, a, b)
    @inbounds for i in eachindex(y, a, b)
        y[i] = exp(a[i] * b[i]) * sin(a[i])
    end
end

# 2. @turbo 优化写法
function heavy_calc_turbo!(y, a, b)
    @turbo for i in eachindex(y, a, b)
        y[i] = exp(a[i] * b[i]) * sin(a[i])
    end
end

println("=== 原生 Base 耗时 ===")
@time heavy_calc_base!(y, a, b)

println("\n=== @turbo 耗时 ===")
@time heavy_calc_turbo!(y, a, b)

