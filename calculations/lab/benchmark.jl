using LoopVectorization
using BenchmarkTools

function benchmark_square_lattice()
    L = 300
    N = L * L
    # === 构造真实的 2D 正方晶格邻居表 ===
    # 邻居是上下左右，具有规律性
    neighbors = Matrix{Int}(undef, 4, N)
    for y in 1:L, x in 1:L
        i = (y-1)*L + x
        # 周期性边界条件
        xm = (x==1) ? L : x-1
        xp = (x==L) ? 1 : x+1
        ym = (y==1) ? L : y-1
        yp = (y==L) ? 1 : y+1
        
        neighbors[1, i] = (y-1)*L + xm
        neighbors[2, i] = (y-1)*L + xp
        neighbors[3, i] = (ym-1)*L + x
        neighbors[4, i] = (yp-1)*L + x
    end
    
    # 构造 Vector{Vector} 形式
    adj_list = [neighbors[:, i] for i in 1:N]
    
    S = rand(3, N)
    H_normal = zeros(3, N)
    H_turbo = zeros(3, N)

    # --- 函数定义 (和之前一样) ---
    function calc_normal!(H, S, adj_list)
        @inbounds for i in eachindex(adj_list)
            nbs = adj_list[i]
            hx, hy, hz = 0.0, 0.0, 0.0
            for n in nbs
                hx += S[1, n]
                hy += S[2, n]
                hz += S[3, n]
            end
            H[1, i] = hx; H[2, i] = hy; H[3, i] = hz
        end
    end

    function calc_turbo!(H, S, neighbors)
        @turbo for i in axes(neighbors, 2)
            hx, hy, hz = 0.0, 0.0, 0.0
            n1=neighbors[1,i]; n2=neighbors[2,i]; n3=neighbors[3,i]; n4=neighbors[4,i]
            
            # 手动展开，避免循环开销
            hx = S[1,n1] + S[1,n2] + S[1,n3] + S[1,n4]
            hy = S[2,n1] + S[2,n2] + S[2,n3] + S[2,n4]
            hz = S[3,n1] + S[3,n2] + S[3,n3] + S[3,n4]
            
            H[1, i] = hx; H[2, i] = hy; H[3, i] = hz
        end
    end

    println("=== 2D Square Lattice Benchmark (N=$N) ===")
    print("Normal: ")
    @btime calc_normal!($H_normal, $S, $adj_list)
    print("Turbo:  ")
    @btime calc_turbo!($H_turbo, $S, $neighbors)
end

benchmark_square_lattice()