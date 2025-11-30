# using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../src/src.jl")
include("model.jl")


dataname = "Heisenberg_XCPeriTria/data"

Lx = 40
Ly = 40
Latt = XCPeriTria(Lx,Ly)
@save "$(dataname)/Latt_$(Lx)x$(Ly)" Latt
params = (J = 1.0,)
H = Dict(
    "J1" => (params.J * diagm([1,1,1]), Tuple(neighbor(Latt))),
    "h" => (-[0.0,0.0,0.0],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[sqrt(3)/2,1/2,0.0],((1,),))
)

Obs = Dict(
    "S" => Tuple(tuple.(1:length(Latt))),
    "SS" => (),
)

thmalgo = Dict(
    "N" => 10,
    "W" => 30,
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
W = 2

lsT = Float64[]
let Ttmp = T0
for i in 1:Nt
    to = TimerOutput()
    if Ttmp < T
        println("T = $(Ttmp) < Tf = $(T) at $(i)th SA")
        break
    end
    Ttmp = α * Ttmp
    algo = Dict(
        "β" => 1 / Ttmp,
        "update" => "Heat bath",
    )
    @timeit to "thermalize" thermalize!(ψ,H,algo,thmalgo)
    data = Dict(
        "N" => N,
        "W" => W,
        "T" => Ttmp,
        "E" => Float64[],
        "M" => Float64[],
        "S" => [zeros(3) for _ in 1:length(Latt)]
    )
    @timeit to "calculate" begin
        for j in 1:N 
            for k in 1:W
                update!(ψ,H,algo)
                data["S"] .+= deepcopy(ψ["S"]) / N / W
            end
            push!(data["E"], calculate_E(ψ,H))
            push!(data["M"], norm(sum(ψ["S"])))        
        end     
    end
    @timeit to "manual GC" GC.gc()
    push!(lsT,Ttmp)
    @save "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(round(Ttmp;digits = 8)).jld2" data
    show(to;title = "$(i)/$(ceil(Int64,log(α,T/T0))) - T = $(round(Ttmp;digits = 4))")
    print("\n")
end
end

@save "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT


