include("../src/src.jl")
include("model.jl")


function SD!(ψ::Dict,H::Dict)
    Nc = 100000
    for count in 1:Nc
        # println("SD - $(count)")
        for i in shuffle(1:length(Latt))
            S′ = calc_heff(ψ,H,i) |> x -> x/norm(x)
            δs[i] = norm(ψ["S"][i] - S′)
            ψ["S"][i] = S′
        end
        if sum(δs) < 1e-16
            println("SD Finished at $(count) with Σ|Sᵢ′ - Sᵢ| = $(sum(δs))")
            return ψ
        end
    end
    println("SD Not Finished after $(Nc) with Σ|Sᵢ′ - Sᵢ| = $(sum(δs))")

    return ψ
end

dataname = "Heisenberg_PeriSqua/data/mt"


Lx = 10
Ly = 10
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 2
T = 0.01
params = (J = 1.0,)
H = Dict(
    "J1" => (params.J * diagm([1,1,1]),
         Tuple(neighbor(Latt)),
         Dict(i => collect(map(x -> (x[1] == i ? x[2] : x[1]), filter(x -> i in x,Tuple(neighbor(Latt)))))  for i in 1:length(Latt))),
    "h" => (-[0.0,0.0,0.0],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0,0,1],((1,),))
)

α = 0.9
# Tf = T0*(α)^(ceil(Int64,log(0.9,0.01/2)))
Nt = 1000
N = 0
W = 2
@load "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT
Tf = lsT[end]
@load "$(dataname)/ψ_$(Lx)x$(Ly)_$(params)_$(N)_$(W)_$(round(Tf;digits = 8)).jld2" ψ

hs = Vector{Vector}(undef,length(Latt))
δs = ones(length(Latt))


SD!(ψ,H)
