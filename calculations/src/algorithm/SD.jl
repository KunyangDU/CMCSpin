
function SD!(ψ::SimpleState{D,L},H::Hamiltonian,algo::SDAlgo = SDAlgo();
    showperstep = 500) where {D,L}
    println("SD begin")
    to = TimerOutput()
    for count in 1:algo.N
        @timeit to "update" δs = update!(ψ,H,algo)
        ϵ = mean(δs)
        if ϵ < algo.tol
            println("SD Finished at $(count) with ⟨|Sᵢ′ - Sᵢ|⟩ = $(sum(δs)/L)")
            return ψ
        end
        if mod(count,showperstep) == 0
            show(to;title = "SD - $(count)/$(algo.N)")
            println("\nϵ = $(ϵ), tol = $(algo.tol)")
        end
    end
    println("SD Not Finished after $(algo.N)")

    return ψ
end

function update!(ψ::SimpleState{D,L},H::Hamiltonian,::SDAlgo) where {D,L}
    δs = zeros(L)
    Nthr = get_num_threads_julia()
    if Nthr > 1 && !isnothing(H.group)
        for sites in H.group
            L′ = length(sites)
            counter = Threads.Atomic{Int64}(1)
            Threads.@sync for _ in 1:Nthr
                Threads.@spawn while true
                    ct = Threads.atomic_add!(counter, 1)
                    ct > L′ && break
                    δs[sites[ct]] = _SD_update_work!(ψ,H,sites[ct])
                end
            end
        end
    else
        for i in shuffle(1:L)
            δs[i] = _SD_update_work!(ψ,H,i)
        end
    end
    return δs
end

function _SD_update_work!(ψ::SimpleState{D,L}, H::Hamiltonian, i::Int64) where {D,L}
    S′ = gradient(ψ,H,i) |> x -> x/norm(x)
    δ = norm(ψ[i] - S′)
    ψ[i] = S′
    return δ
end