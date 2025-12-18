include("../src/src.jl")


Lx = 2
Ly = 2
Latt = PeriSqua(Lx,Ly)

H = Hamiltonian()
addIntr2!(H,Latt,"J1",-diagm([1,1,1]),
# (((1,1),[1,0]),((1,1),[-1,0]),((1,1),[0,1]),((1,1),[0,-1]))
Tuple(neighbor_pbc(Latt))
)
# setgroup!(H,Latt.group)

ψ = SimpleState(hcat([[0.0,1.0,1.0]/sqrt(2) * (1) ^ sum(Latt[i][2]) for i in 1:length(Latt)]...))

vpath,rpath,rnode = vrange([[0,0],[pi,0],[pi,pi],[0,0]],121)
band,ΔE,weight = LSW(ψ,H,vpath;isweight = true)
mask = weight .> 2
# band = band .* mask
fig = Figure()
ax = Axis(fig[1,1],
xticks = (rnode,[L"\Gamma",L"X",L"M",L"\Gamma"]))


for i in 1:length(ψ)
    lines!(ax, rpath, band[i,:],color = :grey)
    # scatter!(ax, rpath, band[i,:],color = get(colorschemes[:Reds],weight[i,:],(0.0,maximum(weight))))
    # scatter!(ax, rpath, band[i,:],markersize = 1 * weight[i,:])
end 

display(fig)
weight'
# band
