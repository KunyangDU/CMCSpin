using LatticeUtilities,LinearAlgebra,CairoMakie,LaTeXStrings,JLD2
include("../lattice/lattice.jl")

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


Lx = 60
Ly = 1
Latt = PeriSqua(Lx,Ly)
params = (J = 1.0,)

nb = filter(x -> x[1] ≠ x[2], neighbor(Latt))
H = Dict(
    "J1" => (params.J * diagm(ones(3)), Tuple(nb)),
    # "h" => (-[hx,hy,hz],Tuple(tuple.(1:length(Latt)))),
    # "hsb" => (-[0.0,0.0,0.0],((1,),))
)


Nt = 600
τ = 0.05
lst = range(0,Nt*τ,Nt+1)


lsq = 2pi*(0:div(Lx,2))/(Lx+1)
lsω = range(0, 2.5, div(Lx,1))
# ω = 0.0

dyS = let dyS = zeros(length(lsq),length(lsω))


for i in 1:20
    tmpdyS = zeros(length(lsq),length(lsω))
    ψ = Dict(
        "S" => [[0,0,(-1)^sum(Latt[i][2])] + 0.03*randn(3) |> x -> x/norm(x) for i in 1:length(Latt)]
    )

    lsS = [deepcopy(ψ["S"]),]
    let S = deepcopy(lsS[1])
        for i in 1:Nt 
            S′ = deepcopy(S)
            for j in eachindex(S)

                heff = zeros(3)

                Jij,sites = H["J1"]
                sites = filter(x -> j ∈ x, (sites))
                nbsites = collect(map(x -> j == x[1] ? x[2] : x[1], sites))
                heff += -sum([Jij * S[s] for s in nbsites])

                S′[j] = rk4(S[j], τ, x -> -cross(S[j],heff)) |> x -> x/norm(x)
            end
            push!(lsS, S′)
            S = S′
        end
    end

    for (iq,q) in enumerate(lsq)
        eq = [exp(1im * dot([q,0], coordinate(Latt,i))) for i in 1:length(Latt)] / sqrt(length(Latt))
        Smq0 = sum(eq' .* lsS[1])
        for (iω,ω) in enumerate(lsω)
            lsSq = map( x -> sum(eq .* x), exp.(1im * ω * lst) .* window.(lst/lst[end]) .* lsS) * τ
            tmpdyS[iq,iω] = abs(sum(map(x -> dot(x,Smq0), lsSq)))
        end
    end

    dyS += tmpdyS
end
dyS
end



@save "dynamics/data/dyS_$(Lx)x$(Ly)_$(Nt)_$(τ)_$(params)_examples.jld2" dyS
