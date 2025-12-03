
function SD!(ψ::SimpleState{D,L},H::Hamiltonian,algo::SDAlgo = SDAlgo();
    showperstep = 500) where {D,L}
    println("SD begin")
    δs = zeros(L)
    to = TimerOutput()
    for count in 1:algo.N
        @timeit to "update" begin 
            for i in shuffle(1:L)
                S′ = gradient(ψ,H,i) |> x -> x/norm(x)
                δs[i] = norm(ψ[i] - S′)
                ψ[i] = S′
            end
        end
        ϵ = sum(δs)/L
        if ϵ < algo.tol
            println("SD Finished at $(count) with ⟨|Sᵢ′ - Sᵢ|⟩ = $(sum(δs)/L)")
            return ψ
        end
        if mod(count,showperstep) == 0
            show(to;title = "SD - $(count)/$(algo.N)")
            println("\nϵ = $(ϵ), tol = $(algo.tol)")
        end
    end
    println("SD Not Finished after $(algo.N) with Σ|Sᵢ′ - Sᵢ| = $(sum(δs))")

    return ψ
end