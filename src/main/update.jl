
function updata_S_Matropolis(S::Vector,η::Number)
    R = [rand() - 0.5 for i in 1:3]
    return S + η*R |> x -> x/norm(x)
end

# function update_Matropolis!(ψ::Dict, H::Dict, algo::Dict)
#     for i in shuffle(1:length(Latt))
        
#         heff = zeros(3)

#         Jij,sites = H["J1"]
#         sites = filter(x -> i ∈ x, (sites))
#         nbsites = collect(map(x -> i == x[1] ? x[2] : x[1], sites))
#         heff += -sum([Jij * ψ["S"][s] for s in nbsites])

#         heff += -H["h"][1]

#         if i in H["hsb"][2][1]
#             heff += -H["hsb"][1]
#         end

#         S = ψ["S"][i]
#         S′ = updata_S(S,algo["η"])

#         ΔE = -dot((S′-S),heff)

#         rand() < min(exp(-algo["β"] * ΔE),1) && (ψ["S"][i] = S′)
#     end
# end

function updata_S_Heatbath(h::Vector,algo::Dict)
    b = algo["β"] * norm(h)

    ŵ = h / norm(h)
    if abs(ŵ[1]) < 0.6 && abs(ŵ[2]) < 0.6
        û = [0.0, -ŵ[3], ŵ[2]]
    else
        û = [-ŵ[2], ŵ[1], 0.0]
    end
    û /= norm(û)
    v̂ = cross(ŵ,û)

    cθ = max(min(log(exp(-b) + rand()*(exp(b) - exp(-b))) / b,1),-1)
    sθ = sqrt(1 - cθ^2)
    sϕ,cϕ = sincos(2pi * rand())

    return hcat(û,v̂,ŵ) * [sθ*cϕ,sθ*sϕ,cθ]
    # return hcat(ŵ,nullspace(ŵ')) * [sθ*cϕ,sθ*sϕ,cθ]
end

function update!(ψ::Dict, H::Dict, algo::Dict)
    # Jij,_,nbsites = H["J1"]
    for i in shuffle(1:length(Latt))
        heff = calc_heff(ψ,H,i)

        if algo["update"] == "Heat bath"
            ψ["S"][i] = updata_S_Heatbath(heff,algo)
        elseif algo["update"] == "Matropolis"
            S = ψ["S"][i]
            S′ = updata_S_Matropolis(S,algo["η"])
            ΔE = -dot((S′-S),heff)
            rand() < min(exp(-algo["β"] * ΔE),1) && (ψ["S"][i] = S′)
        end
    end
end

function calc_heff(ψ::Dict, H::Dict, i::Int64)
    heff = zeros(3)
    heff += - H["J1"][1] * sum(ψ["S"][H["J1"][3][i]])
    heff += -H["h"][1]
    if i in H["hsb"][2][1]
        heff += -H["hsb"][1]
    end
    return heff
end

# function update!(ψ::Dict, H::Dict, algo::Dict)
#     tasklist = shuffle(1:length(Latt))
#     Nthr = get_num_threads_julia()
#     if Nthr > 1
#         Lock = Threads.ReentrantLock()
#         counter = Threads.Atomic{Int64}(1)
#         Threads.@sync for _ in 1:Nthr
#             Threads.@spawn while true
#                 ct = Threads.atomic_add!(counter, 1)
#                 ct > length(tasklist) && break
#                 i = tasklist[ct]

#                 heff = zeros(3)
#                 heff += -H["h"][1]
#                 if i in H["hsb"][2][1]
#                     heff += -H["hsb"][1]
#                 end

#                 Jij,sites = H["J1"]
#                 sites = filter(x -> i ∈ x, (sites))
#                 nbsites = collect(map(x -> i == x[1] ? x[2] : x[1], sites))

#                 lock(Lock)
#                 try
#                     heff += -sum([Jij * ψ["S"][s] for s in nbsites])
#                     if algo["update"] == "Heat bath"
#                         ψ["S"][i] = updata_S_Heatbath(heff,algo)
#                     elseif algo["update"] == "Matropolis"
#                         S = ψ["S"][i]
#                         S′ = updata_S_Matropolis(S,algo["η"])
#                         ΔE = -dot((S′-S),heff)
#                         rand() < min(exp(-algo["β"] * ΔE),1) && (ψ["S"][i] = S′)
#                     end
#                 catch
#                     rethrow()
#                 finally
#                     unlock(Lock)
#                 end
#             end
#         end
#     else
#         for i in shuffle(1:length(Latt))
#             heff = zeros(3)

#             Jij,sites = H["J1"]
#             sites = filter(x -> i ∈ x, (sites))
#             nbsites = collect(map(x -> i == x[1] ? x[2] : x[1], sites))
#             heff += -sum([Jij * ψ["S"][s] for s in nbsites])

#             heff += -H["h"][1]

#             if i in H["hsb"][2][1]
#                 heff += -H["hsb"][1]
#             end

#             if algo["update"] == "Heat bath"
#                 ψ["S"][i] = updata_S_Heatbath(heff,algo)
#             elseif algo["update"] == "Matropolis"
#                 S = ψ["S"][i]
#                 S′ = updata_S_Matropolis(S,algo["η"])
#                 ΔE = -dot((S′-S),heff)
#                 rand() < min(exp(-algo["β"] * ΔE),1) && (ψ["S"][i] = S′)
#             end
#         end
#     end
# end