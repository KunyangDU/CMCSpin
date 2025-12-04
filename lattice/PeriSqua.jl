using Random
include("../src/src.jl")


Lx = 10
Ly = 4

Latt = PeriSqua(Lx,Ly,Tuple(shuffle(1:length(ψ))))

figsize = (height = (Ly+1)*50, width = (Lx + 1)*50)

fig = Figure()
ax = Axis(fig[1,1];autolimitaspect = true,figsize...)

# latticescatter!(ax,Latt)
plotLatt!(ax,Latt;site = true,tplevel = (1,))
resize_to_layout!(fig)
display(fig)

neighbor(Latt,1)

save("lattice/figures/PeriSqua_$(Lx)x$(Ly).pdf",fig)
save("lattice/figures/PeriSqua_$(Lx)x$(Ly).png",fig)

# Latt.unitcell.lattice_vecs * neighbor_pbc(Latt)[10][2]
neighbor_pbc(Latt)
