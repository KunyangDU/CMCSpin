include("../src/src.jl")

Lx = 10
Ly = 10

Latt = PeriSqua(Lx,Ly)

groups = group(Latt,[1,],2)
Latt = PeriSqua(Lx,Ly,reverse_map(Tuple(vcat(groups...))))
groups = map(g -> Tuple(Latt.sitemap[g]),groups)
Latt.group = Tuple(groups)
@save "lab/data/Latt_$(Lx)x$(Ly).jld2" Latt

T0 = 2.0
Tf = 0.5
α = 0.9
params = (
    J = 1.0,
)

θ = 0.0 * pi
ϕ = 0.0 * pi


H = Hamiltonian()
addIntr2!(H,Latt,"J1",params.J * diagm([1,1,1]),Tuple(neighbor(Latt)))
addIntr1!(H,Latt,"H",-[0,0,0],Tuple(1:length(Latt)))
# addIntr1!(H,Latt,"h",-[0.0,0.0,0.0],(1,))
setgroup!(H,Latt.group)

lsHf = 0:0.1:8.0

for Hf in lsHf

Hx = Hf*sin(θ)*cos(ϕ)
Hy = Hf*sin(θ)*sin(ϕ)
Hz = Hf*cos(θ) 

H.param[2] = -[Hx,Hy,Hz]

ψ = SimpleState(Latt)
normalize!(ψ)
_,lsT,lsE,_ = SA!(ψ,H,SAAlgo(T0,Tf,α,ThmAlgo(50,1000,1.0)))
SD!(ψ,H)

@save "lab/data/ψ_$(Lx)x$(Ly)_$(round(θ/pi))_$(round(ϕ/pi))_$(Hf)_$(params)_$(T0)_$(Tf)_$(α).jld2" ψ
end