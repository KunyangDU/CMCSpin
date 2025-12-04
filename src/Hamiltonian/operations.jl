
function setgroup!(H::Hamiltonian, group::Tuple)
    H.group = group
    for rel in H.relative
        isnothing(rel) && continue
        for g in H.group,i in g 
            @assert !issubset(rel[i], g) "group divided wrong"
        end
    end
    return H
end

# function measure(ψ::SimpleState, H::Hamiltonian)
#     E = 0.0
#     Nthr = get_num_threads_julia()
#     for i in eachindex(H.name)
#         L = length(H.site[i])
#         if Nthr > 1
#             Es = zeros(Nthr)
#             Threads.@sync for iNthr in 1:Nthr
#                 Threads.@spawn for j in get_chunk_range(L,Nthr,iNthr)
#                     Es[Threads.threadid()] += _micro_measure!(ψ,H.param[i],H.site[i][j])
#                 end
#             end
#             E += sum(Es)
#         else
#             for s in H.site[i]
#                 E += _micro_measure!(ψ,H.param[i],s)
#             end
#         end
#     end
#     return E
# end

function measure(ψ::SimpleState, H::Hamiltonian)
    E = 0.0
    Nthr = get_num_threads_julia()
    for i in eachindex(H.name)
        L = length(H.site[i])
        if Nthr > 1
            Es = zeros(Nthr)
            counter = Threads.Atomic{Int64}(1)
            Threads.@sync for _ in 1:Nthr
                Threads.@spawn while true
                    ct = Threads.atomic_add!(counter, 1)
                    ct > L && break
                    Es[Threads.threadid()] += _micro_measure!(ψ,H.param[i],H.site[i][ct])
                end
            end
            E += sum(Es)
        else
            for s in H.site[i]
                E += _micro_measure!(ψ,H.param[i],s)
            end
        end
    end
    return E
end

function _micro_measure!(ψ::SimpleState,P::Vector,i::Int64)
    return P' * ψ[i]
end

function _micro_measure!(ψ::SimpleState,P::Matrix,i::NTuple{2,Int64})
    return ψ[i[1]]' * P * ψ[i[2]]
end
