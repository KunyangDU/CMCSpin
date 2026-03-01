
using Clustering,MultivariateStats
include("../src/src.jl")


Lx = 10
Ly = 10
@load "lab/data/Latt_$(Lx)x$(Ly).jld2" Latt

T0 = 2.0
Tf = 0.5
α = 0.9
params = (
    J = 1.0,
)

θ = 0.0 * pi
ϕ = 0.0 * pi
lsHf = 0.0:0.2:8.0
SSF = zeros(length(Latt),length(lsHf))
Hf = 0.0
@load "lab/data/ψ_$(Lx)x$(Ly)_$(round(θ/pi))_$(round(ϕ/pi))_$(Hf)_$(params)_$(T0)_$(Tf)_$(α).jld2" ψ

Latt′,ψ′ = minimum_pattern(Latt,ψ)

ψ′.pattern


