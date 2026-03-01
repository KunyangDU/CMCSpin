include("../../src/src.jl")

dataname = "examples/KitaevGamma"

Lx = 12
Ly = 12

params = (J = 1.0,K = 0.0,Γ = 0.0,Γ′ = 0.0)

θ = 0.0 * pi
ϕ = 0.0 * pi

T0 = 2.0
Tf = 0.01
α = 0.95

for Hf in 0.0:1.0:8.0



@load "$(dataname)/data/Latt_$(Lx)x$(Ly).jld2" Latt
@load "$(dataname)/data/ψ_$(params)_$(θ)_$(ϕ)_$(Hf)_$(T0)_$(Tf)_$(α).jld2" ψ


data = Dict(
    "S" => ψ.pattern
)

@save "$(dataname)/data/data_$(params)_$(θ)_$(ϕ)_$(Hf)_$(T0)_$(Tf)_$(α).jld2" data

end