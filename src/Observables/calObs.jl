function calculate_E(ψ::Dict, H::Dict)
    E = 0.0

    Jij,sites,_ = H["J1"]
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