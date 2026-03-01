include("../src/src.jl")

Lx = 10
Ly = 10

Latt = PeriSqua(Lx,Ly)
groups = group(Latt,[1,],2)
Latt = PeriSqua(Lx,Ly,reverse_map(Tuple(vcat(groups...))))
groups = map(g -> Tuple(Latt.sitemap[g]),groups)
Latt.group = Tuple(groups)

θ = 0.0 * pi
ϕ = 0.0 * pi
Hf = 1.0

T0 = 2.0
Tf = 0.5
α = 0.9

Hx,Hy,Hz = Hf*sin(θ)*cos(ϕ), Hf*sin(θ)*sin(ϕ), Hf*cos(θ) 

H = Hamiltonian()
addIntr2!(H,Latt,"J1",diagm([1,1,1]),Tuple(neighbor(Latt)))
addIntr1!(H,Latt,"H",-[Hx,Hy,Hz],Tuple(1:length(Latt)))
addIntr1!(H,Latt,"h",-[0.0,0.0,1.0],(1,))
setgroup!(H,Latt.group)

ψ = SimpleState(Latt)
normalize!(ψ)
_,lsT,lsE,_ = SA!(ψ,H,SAAlgo(T0,Tf,α,ThmAlgo(50,1000,1.0)))
SD!(ψ,H)


