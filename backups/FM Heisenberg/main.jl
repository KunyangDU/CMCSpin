using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")

dataname = "FM Heisenberg/data"

Lx = 10
Ly = 10
Latt = PeriSqua(Lx,Ly)
@save "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
# η = 0.02
N = 100000
n = 100


H = Dict(
    "J1" => (-diagm(ones(3)), Tuple(neighbor(Latt))),
    "h" => (-[hx,hy,hz],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0.0,0.0,0.0],((1,),))
)

Obs = Dict(
    "S" => Tuple(tuple.(1:length(Latt))),
    "SS" => (),
    # "SS" => Tuple([(i,j) for i in 1:length(Latt) for j in i+1:length(Latt)])
)

for T in 1.0:0.1:2.0
algo = Dict(
    "β" => 1 / T,
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

@save "$(dataname)/total_data_Heatbath_$(Lx)x$(Ly)_$(T)_$(N)_$(n)" total_data
end


