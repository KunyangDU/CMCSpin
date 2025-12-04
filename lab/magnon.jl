include("../src/src.jl")


Lx = 2
Ly = 2
Latt = PeriSqua(Lx,Ly)

H = Hamiltonian()
addIntr2!(H,Latt,"J1",diagm([1,1,1]),Tuple(
neighbor_pbc(Latt)))
# setgroup!(H,Latt.group)

ψ = SimpleState(hcat([[0.0,0.0,1.0] * (-1) ^ sum(Latt[i][2]) for i in 1:length(Latt)]...))

vpath,rpath,rnode = vrange([[0,0],[0,pi],[pi,pi],[0,0]],127)
band,ΔE = LSW(ψ,H,vpath;isweight = false)
# ΔE = sum(band,dims = 1)[:] - As

fig = Figure()
ax = Axis(fig[1,1];
width = 300, height = 200,
ylabel = L"\hbar \omega",
xticks = (rnode,[L"(0,0)",L"(0,\pi)",L"(\pi,\pi)",L"(0,0)"]))
for i in 1:length(ψ)
    lines!(ax, rpath, band[i,:],color = :black)
    # scatter!(ax, rpath, band[i,:],color = get(colorschemes[:Reds],weight[i,:],(0.0,maximum(weight))))
    # scatter!(ax, rpath, band[i,:],markersize = 1 * weight[i,:])
end 

xlims!(ax, extrema(rpath))
ylims!(ax, 0, maximum(band) * 1.2)
resize_to_layout!(fig)
display(fig)

ΔE