
include("../src/src.jl")


Lx = 6
Ly = 6

Latt = PeriSqua(Lx,Ly)

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
