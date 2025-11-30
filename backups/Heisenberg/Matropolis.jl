using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")



# function ()
    
# end


Lx = 6
Ly = 6
Latt = PeriSqua(Lx,Ly)
@save "Heisenberg/data/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
η = 0.02
N = 10000
n = 100


H = Dict(
    "J1" => (diagm(ones(3)), Tuple(neighbor(Latt))),
    "h" => (-[hx,hy,hz],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0.0,0.0,1],((1,),))
)

Obs = Dict(
    "S" => Tuple(tuple.(1:length(Latt))),
    "SS" => Tuple([(i,j) for i in 1:length(Latt) for j in i+1:length(Latt)])
)

for β in 10

algo = Dict(
    "β" => β,
    "update" => "Matropolis",
    "η" => η 
)

ψ = Dict(
    "S" => [randn(3) |> x -> x/norm(x) for i in 1:length(Latt)]
)


# total_data = Vector{Any}(undef,div(N,n))
total_data = []
lsE = zeros(n)
lsdata = Vector(undef,n)
@time for i in 1:N
    update!(ψ,H,algo)
    lsE[mod(i - 1, n) + 1] = calculate_E(ψ,H)
    # lsdata[mod(i - 1, n) + 1] = calObs(ψ,Obs)
    if mod(i,n) == 0
        data = calObs(ψ,Obs)
        data["E"] = mean(lsE)
        data["σE"] = std(lsE)
        # data["SS"] = let 
        #     tmpdata = Dict()
        #     for key in keys(lsdata[1]["SS"])
        #         tmpdata[key] = 0.0
        #         for i in 1:n
        #             tmpdata[key] += lsdata[i]["SS"][key] / n
        #         end
        #     end
        #     tmpdata
        # end
        # data["S"] = let 
        #     tmpdata = Dict()
        #     for key in keys(lsdata[1]["S"])
        #         tmpdata[key] = zeros(3)
        #         for i in 1:n
        #             tmpdata[key] += lsdata[i]["S"][key] / n
        #         end
        #     end
        #     tmpdata
        # end
        push!(total_data, data)
    end
end

@save "Heisenberg/data/total_data_Matropolis_$(Lx)x$(Ly)_$(β)_$(η)_$(N)_$(n)" total_data
end


