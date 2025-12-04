include("../src/src.jl")

Lx = 40
Ly = 40

Latt = PeriSqua(Lx,Ly)
groups = group(Latt,[1,],2)
Latt = PeriSqua(Lx,Ly,reverse_map(Tuple(vcat(groups...))))
groups = map(g -> Tuple(Latt.sitemap[g]),groups)
Latt.group = Tuple(groups)

H = Hamiltonian()
addIntr2!(H,Latt,"J1",-diagm([1,1,1]),Tuple(neighbor(Latt)))
addIntr1!(H,Latt,"H",-[0.0,0.0,0.1],Tuple(1:length(Latt)))
addIntr1!(H,Latt,"h",-[0.0,0.0,1.0],(1,))
setgroup!(H,Latt.group)

ψ = SimpleState(Latt)
normalize!(ψ)
T0 = 2.0
Tf = 0.5
α = 0.9
_,lsT,lsE,_ = SA!(ψ,H,SAAlgo(T0,Tf,α,ThmAlgo(50,1000,1.0)))
SD!(ψ,H)
# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax,lsE)
# display(fig)
@show measure(ψ,H)
ψ.pattern

