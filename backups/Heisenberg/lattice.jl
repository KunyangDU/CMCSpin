using JLD2,LatticeUtilities,CairoMakie
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")


Lx = 6
Ly = 6

Latt = PeriSqua(Lx,Ly)

figsize = (height = (Ly+1)*50, width = (Lx + 1)*50)

fig = Figure()
ax = Axis(fig[1,1];autolimitaspect = true,figsize...)

# latticescatter!(ax,Latt)
plotLatt!(ax,Latt,[1,1];site = true)
resize_to_layout!(fig)
display(fig)

neighbor(Latt,1)