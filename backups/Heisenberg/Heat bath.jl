using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")

Lx = 12
Ly = 12
Latt = PeriSqua(Lx,Ly)
@save "Heisenberg/data/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
# η = 0.02
N = 50000
n = 200


H = Dict(
    "J1" => (diagm(ones(3)), Tuple(neighbor(Latt))),
    "h" => (-[hx,hy,hz],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0.0,0.0,1],((1,),))
)

Obs = Dict(
    "S" => Tuple(tuple.(1:length(Latt))),
    "SS" => Tuple([(i,j) for i in 1:length(Latt) for j in i+1:length(Latt)])
)

for β in 1:1.0:20.0
algo = Dict(
    "β" => β,
    "update" => "Heat bath"
)
ψ = Dict(
    "S" => [randn(3) |> x -> x/norm(x) for i in 1:length(Latt)]
)


# total_data = Vector{Any}(undef,div(N,n))
total_data = []
# lsE = zeros(n)
# lsdata = Vector(undef,n)
@time for i in 1:N
    update!(ψ,H,algo)
    if mod(i,n) == 0
        data = calObs(ψ,Obs)
        data["E"] = calculate_E(ψ,H)
        push!(total_data, data)
    end
end

@save "Heisenberg/data/total_data_Heatbath_$(Lx)x$(Ly)_$(β)_$(N)_$(n)" total_data
end


