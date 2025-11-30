using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")


dataname = "analysis/data"

Lx = 8
Ly = 8
Latt = PeriSqua(Lx,Ly)
@save "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

H = Dict(
    "J1" => (diagm([1,1,1]), Tuple(neighbor(Latt))),
    "h" => (-[0.0,0.0,0.0],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0.0,0.0,1],((1,),))
)

Obs = Dict(
    "S" => Tuple(tuple.(1:length(Latt))),
    "SS" => (),
)

thmalgo = Dict(
    "N" => 10,
    "W" => 50,
    "ϵ" => 0.5
)

# BSAalgo = Dict(
#     "N" => 2^15,
#     "W" => 2 .^ (0:1:9)
# )

ψ = Dict(
    "S" => [randn(3) |> x -> x/norm(x) for i in 1:length(Latt)]
)

T0 = 2
T = 0.01

α = 0.9
Nt = 1000
N = 1000

lsT,lsE,lsS = let Ttmp = T0
lsS = Vector[]
lsE = Vector[]
lsT = Float64[]
for i in 1:Nt
    if Ttmp < T
        println("T = $(Ttmp) < Tf = $(T) at $(i)th SA")
        break
    end
    Ttmp = α * Ttmp
    algo = Dict(
        "β" => 1 / Ttmp,
        "update" => "Heat bath",
    )
    _,_,_,to = thermalize!(ψ,H,algo,thmalgo)
    show(to)
    lsEtmp = Vector{Float64}(undef,N)
    lsStmp = Vector{Vector}(undef,N)
    for j in 1:N 
        update!(ψ,H,algo)
        lsEtmp[j] = calculate_E(ψ,H)
        lsStmp[j] = deepcopy(ψ["S"])
    end
    push!(lsE,lsEtmp)
    push!(lsS,lsStmp)
    push!(lsT,Ttmp)
end
lsT,lsE,lsS
end
data = Dict(
    "T" => lsT,
    "E" => lsE,
    "S" => lsS
)

final_data = Dict(
    # "T" => lsT,
    "E" => lsE[end],
    "S" => lsS[end]
)
@save "$(dataname)/dara_SA_$(Lx)x$(Ly)_$(T0)_$(T)_$(α)_$(Nt)_$(N).jld2" data
@save "$(dataname)/final_data_SA_$(Lx)x$(Ly)_$(T0)_$(T)_$(α)_$(Nt)_$(N).jld2" final_data
@save "$(dataname)/ψ_SA_$(Lx)x$(Ly)_$(T0)_$(T)_$(α)_$(Nt)_$(N).jld2" ψ

# N = 1000
# lsE = zeros(N)
# for i in 1:N
#     update!(ψ,H,algo)
#     lsE[i] = calculate_E(ψ,H)
# end
# mean(lsE) / length(Latt)
# lines(lsE)
# lines(lsT,lsE)

# fig = Figure()
# ax = Axis(fig[1,1],xscale = log10)

# scatterlines!(ax,lsT,lsE / length(Latt))
# display(fig)


# 1