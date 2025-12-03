

function Matropolis!(S::Vector,η::Number)
    R = [rand() - 0.5 for i in 1:3]
    return S + η*R |> x -> x/norm(x)
end

function Heatbath!(h::Vector,β::Number,rθ::Number = rand(),rϕ::Number = rand())
    b = β * norm(h)

    ŵ = h / norm(h)
    if abs(ŵ[1]) < 0.6 && abs(ŵ[2]) < 0.6
        û = [0.0, -ŵ[3], ŵ[2]]
    else
        û = [-ŵ[2], ŵ[1], 0.0]
    end
    û /= norm(û)
    v̂ = cross(ŵ,û)

    cθ = max(min(log(exp(-b) + rθ*(exp(b) - exp(-b))) / b,1),-1)
    sθ = sqrt(1 - cθ^2)
    sϕ,cϕ = sincos(2pi * rϕ)

    return hcat(û,v̂,ŵ) * [sθ*cϕ,sθ*sϕ,cθ]
end

function update!(ψ::SimpleState, H::Hamiltonian, algo::UpdateAlgo{Sch}) where Sch
    for i in shuffle(1:length(Latt))
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
    return ψ
end