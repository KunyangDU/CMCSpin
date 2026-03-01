include("../../src/src.jl")
include("model.jl")

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



Latt = XCPH1(Lx,Ly)
@save "$(dataname)/data/Latt_$(Lx)x$(Ly).jld2" Latt

Hv = Hf * (P * [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)] )
H = CubicHamiltonian(Latt;params...,Hv = Hv)

ψ = SimpleState(Latt)
normalize!(ψ)
_,lsT,lsE,_ = SA!(ψ,H,SAAlgo(T0,Tf,α,ThmAlgo(50,1000,1.0)))
SD!(ψ,H,SDAlgo(10000,1e-8))

@save "$(dataname)/data/ψ_$(params)_$(θ)_$(ϕ)_$(Hf)_$(T0)_$(Tf)_$(α).jld2" ψ
end