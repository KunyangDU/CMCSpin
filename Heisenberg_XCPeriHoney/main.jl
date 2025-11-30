# using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../src/src.jl")
include("model.jl")


dataname = "Heisenberg_XCPeriHoney/data"

Lx = 8
Ly = 8
Latt = XCPeriHoney(Lx,Ly)
@save "$(dataname)/Latt_$(Lx)x$(Ly)" Latt
params = (J = -1.0,)
H = Dict(
    "J1" => (params.J * diagm([1,1,1]), Tuple(neighbor(Latt))),
    "h" => (-[0.0,0.0,0.0],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0.0,0.0,1.0],((1,),))
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

ψ = Dict(
    "S" => [randn(3) |> x -> x/norm(x) for i in 1:length(Latt)]
)

T0 = 10
T = 0.01

α = 0.9
Nt = 1000
N = 1000
W = 5

lsT = Float64[]
let Ttmp = T0
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
    thermalize!(ψ,H,algo,thmalgo)
    data = Dict(
        "N" => N,
        "W" => W,
        "T" => Ttmp,
        "E" => Float64[],
        "M" => Float64[],
        "S" => [zeros(3) for _ in 1:length(Latt)]
    )
    for j in 1:N 
        for k in 1:W
            update!(ψ,H,algo)
            data["S"] .+= deepcopy(ψ["S"]) / N / W
        end
        push!(data["E"], calculate_E(ψ,H))
        push!(data["M"], norm(sum(ψ["S"])))        
    end
    push!(lsT,Ttmp)
    # push!(data["data"], data)
    @save "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(round(Ttmp;digits = 8)).jld2" data
end
end

@save "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT


