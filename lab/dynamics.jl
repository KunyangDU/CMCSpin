include("../src/src.jl")


Lx = 20
Ly = 4
Latt = PeriSqua(Lx,Ly)

H = Hamiltonian()
addIntr2!(H,Latt,"J1",-diagm([1,1,1]),Tuple(
neighbor_pbc(Latt)))
# setgroup!(H,Latt.group)

ψ = SimpleState(hcat([[0.0,0.0,1.0] * (-1) ^ sum(Latt[i][2]) for i in 1:length(Latt)]...))

vpath,rpath,rnode = vrange([[0,pi],[pi,pi]],100)
band,weight = LSW(ψ,H,vpath;isweight = true)

fig = Figure()
ax = Axis(fig[1,1])
for i in 1:length(ψ)
    lines!(ax, rpath, band[i,:],color = :grey)
    # scatter!(ax, rpath, band[i,:],color = get(colorschemes[:Reds],weight[i,:],(0.0,maximum(weight))))
    scatter!(ax, rpath, band[i,:],markersize = 1 * weight[i,:])
end 

display(fig)
