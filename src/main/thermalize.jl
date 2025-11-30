function thermalize!(ψ::Dict, H::Dict, algo::Dict, thmalgo::Dict)
    lsE = Float64[]
    lsM = Float64[]
    E = Inf
    M = Inf
    N = thmalgo["N"]
    W = thmalgo["W"]
    to = TimerOutput()

    for i in 1:N
        tmplsE = zeros(W)
        tmplsM = zeros(W)
        for j in 1:W 
            @timeit to "update!" update!(ψ,H,algo)
            tmplsE[j] = calculate_E(ψ,H)
            tmplsM[j] = norm(mean(ψ["S"]))
            push!(lsE,tmplsE[j])
            push!(lsM,tmplsM[j])
        end
        E′ = mean(tmplsE)
        M′ = mean(tmplsM)
        if abs(E - E′) < thmalgo["ϵ"] * std(tmplsE) / sqrt(N/W) && abs(M - M′) < thmalgo["ϵ"] * std(tmplsM) / sqrt(N/W)
            println("THERMALIZED at $(i)th iteration: |Q′ - Q| < $(thmalgo["ϵ"])σ/sqrt(W)")
            return ψ,lsE,lsM,to
        end
        if i == N 
            println("NOT THERMALIZED with $(N) iterations: |E′ - E| = $(abs(E - E′)), σ = $(std(tmplsE))")
        end
        E = E′
        M = M′
    end

    return ψ,lsE,lsM,to
end