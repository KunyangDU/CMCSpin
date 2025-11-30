
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

    cθ = log(exp(-b) + rand()*(exp(b) - exp(-b))) / b
    sθ = sqrt(1 - cθ^2)
    sϕ,cϕ = sincos(2pi * rand())
    S′ = [sθ*cϕ,sθ*sϕ,cθ]

    ŵ = h / norm(h)
    if abs(ŵ[1]) < 0.6 && abs(ŵ[2]) < 0.6
        û = [0.0, -ŵ[3], ŵ[2]]
    else
        û = [-ŵ[2], ŵ[1], 0.0]
    end
    û /= norm(û)
    v̂ = cross(ŵ,û)

    return hcat(û,v̂,ŵ) * S′
end

function update!(ψ::Dict, H::Dict, algo::Dict)
    for i in shuffle(1:length(Latt))
        heff = zeros(3)

        Jij,sites = H["J1"]
        sites = filter(x -> i ∈ x, (sites))
        nbsites = collect(map(x -> i == x[1] ? x[2] : x[1], sites))
        heff += -sum([Jij * ψ["S"][s] for s in nbsites])

        heff += -H["h"][1]

        if i in H["hsb"][2][1]
            heff += -H["hsb"][1]
        end

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

function calculate_E(ψ::Dict, H::Dict)
    E = 0.0

    Jij,sites = H["J1"]
    for (i,j) in sites
        E += ψ["S"][i]' * Jij * ψ["S"][j]
    end

    hi,sites = H["h"]
    for (i,) in sites
        E += hi' * ψ["S"][i]
    end

    hi,sites = H["hsb"]
    for (i,) in sites
        E += hi' * ψ["S"][i]
    end

    return E
end

function calObs(ψ::Dict, Obs::Dict, H::Union{Nothing, Dict} = nothing)
    data = Dict()
    data["SS"] = Dict()
    data["S"] = Dict()
    for (i,j) in Obs["SS"]
        data["SS"][(i,j)] = dot(ψ["S"][i], ψ["S"][j])
    end

    for (i,) in Obs["S"]
        data["S"][(i,)] = ψ["S"][i]
    end

    # if !isnothing(H)
    #     data["E"] = calculate_E(ψ,H)
    # end

    return data
end

function FT2(Latt::SimpleLattice, SM::Matrix, k::Tuple)
    N = length(Latt)
    @assert SM ≈ SM'
    # k = [0.0,0.0]
    SSF = 0.0
    for i in 1:N, j in i+1:N
        SSF += SM[i,j] * exp(1im*dot(collect(k),coordinate(Latt,i) - coordinate(Latt,j))) / N
    end
    for i in 1:N
        SSF += SM[i,i] / N
    end
    return real(SSF) / N
end

function FT2(Latt::SimpleLattice, SM::Matrix, tk::Vector)
    return [FT2(Latt,SM,k) for k in tk]
end

function rk4(y::Union{Number, Vector}, dt::Number, f::Function)
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    
    return @. y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
end

function window(x::Number)
    if 1/2 < x ≤ 1
        return -2*(x-1)^3
    elseif 0 ≤ x ≤ 1/2
        return 1 - 6x^2 + 6x^3
    elseif -1/2 ≤ x < 0
        return 1 - 6x^2 - 6x^3
    elseif -1 ≤ x < -1/2
        return 2*(x+1)^3
    else
        return 0
    end
end


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

function BSA!(ψ::Dict,H::Dict,BSAalgo::Dict)
    N = BSAalgo["N"]
    lsE = Float64[]

    for i in 1:N
        update!(ψ,H,algo)
        push!(lsE, calculate_E(ψ,H))
    end

    lsW = BSAalgo["W"]
    lsK = div.(N,lsW)

    return BSA(lsE,lsW,lsK)
end

function BSA(lsQ::Vector,lsW::Vector,lsK::Vector)
    lsσ = zeros(length(lsW))
    for (i,W) in enumerate(lsW)
        K = div(N,W)
        lsσ[i] = std([mean(lsQ[(j-1)*W + 1:j*W]) for j in 1:K]) / sqrt(K)
    end

    lsδ = lsσ ./ (2*sqrt.(lsK .- 1))
    τ = let A = 0, B = 0
        for i in 1:length(lsσ)-1
            if abs(lsσ[i + 1] - lsσ[i]) < lsδ[i + 1] + lsδ[i]
                A += lsσ[i] / lsδ[i]^2
                B += 1 / lsδ[i]^2
            end
        end
        ((A/B/lsσ[1])^2 - 1)/2
    end

    return τ
end


