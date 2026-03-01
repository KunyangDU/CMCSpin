
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


