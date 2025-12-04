

function Matropolis!(S::Vector,η::Number)
    R = [rand() - 0.5 for i in 1:3]
    return S + η*R |> x -> x/norm(x)
end

function Heatbath!(h::Vector,β::Number,rθ::Number = rand(),rϕ::Number = rand())
    b = β * norm(h)

    û,v̂,ŵ = _local_axis(h)

    cθ = max(min(log(exp(-b) + rθ*(exp(b) - exp(-b))) / b,1),-1)
    sθ = sqrt(1 - cθ^2)
    sϕ,cϕ = sincos(2pi * rϕ)

    return hcat(û,v̂,ŵ) * [sθ*cϕ,sθ*sϕ,cθ]
end

# function update!(ψ::SimpleState, H::Hamiltonian, algo::UpdateAlgo{Sch}) where Sch
#     Nthr = get_num_threads_julia()
#     if Nthr > 1 && !isnothing(H.group)
#         for sites in H.group
#             L = length(sites)
#             Threads.@sync for iNthr in 1:Nthr
#                 Threads.@spawn begin
#                     for i in shuffle(get_chunk_range(L,Nthr,iNthr))
#                         _update_work!(ψ,H,algo,sites[i])
#                     end
#                 end
#             end
#         end
#     else
#         for i in shuffle(1:length(Latt))
#             _update_work!(ψ,H,algo,i)
#         end
#     end
#     return ψ
# end

function update!(ψ::SimpleState, H::Hamiltonian, algo::UpdateAlgo{Sch}) where Sch
    Nthr = get_num_threads_julia()
    if Nthr > 1 && !isnothing(H.group)
        for sites in H.group
            L = length(sites)
            counter = Threads.Atomic{Int64}(1)
            Threads.@sync for _ in 1:Nthr
                Threads.@spawn while true
                    ct = Threads.atomic_add!(counter, 1)
                    ct > L && break
                    _update_work!(ψ,H,algo,sites[ct])
                end
            end
        end
    else
        for i in shuffle(1:length(Latt))
            _update_work!(ψ,H,algo,i)
        end
    end
    return ψ
end

function _update_work!(ψ::SimpleState{D,L}, H::Hamiltonian, algo::UpdateAlgo{Sch}, i::Int64) where {Sch,D,L}
    heff = gradient(ψ,H,i)
    if Sch == HeatBath
        ψ[i] = Heatbath!(heff,1 / algo.T)
    elseif Sch == Matropolis
        S = ψ[i]
        S′ = Matropolis!(S,algo.scheme.η)
        ΔE = -dot((S′-S),heff)
        rand() < min(exp(-ΔE/algo.T),1) && (ψ[i] = S′)
    end
end
