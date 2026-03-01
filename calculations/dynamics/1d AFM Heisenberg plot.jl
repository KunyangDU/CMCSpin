using LatticeUtilities,LinearAlgebra,CairoMakie,LaTeXStrings,JLD2
include("../lattice/lattice.jl")

Lx = 60
Ly = 1
Latt = PeriSqua(Lx,Ly)
params = (J = 1.0,)

nb = filter(x -> x[1] ≠ x[2], neighbor(Latt))
H = Dict(
    "J1" => (params.J * diagm(ones(3)), Tuple(nb)),
    # "h" => (-[hx,hy,hz],Tuple(tuple.(1:length(Latt)))),
    # "hsb" => (-[0.0,0.0,0.0],((1,),))
)


Nt = 600
τ = 0.05
lst = range(0,Nt*τ,Nt+1)


lsq = 2pi*(0:div(Lx,2))/(Lx+1)
lsω = range(0, 2.5, div(Lx,1))
# ω = 0.0

@load "dynamics/data/dyS_$(Lx)x$(Ly)_$(Nt)_$(τ)_$(params)_examples.jld2" dyS


fig = Figure()
ax = Axis(fig[1,1];
title = "$(Lx)x$(Ly) AFM Heisenberg, ϵ=$(round(8/(Nt*τ);digits = 3))",
width = 250,height = 300,
xticks = ((0:0.5:1)*pi,[L"0",L"\pi/2",L"\pi"]),
ylabel = L"\omega",
xlabel = L"q")

hm = heatmap!(ax, lsq,lsω, dyS, colorrange = (0,1))
Colorbar(fig[1,2],hm,label = L"S(q,\omega)")

xlims!(ax,0,pi)
ylims!(ax,extrema(lsω))

resize_to_layout!(fig)

display(fig)


save("dynamics/figures/dynamics_$(Lx)x$(Ly)_$(Nt)_$(τ)_examples.pdf",fig)
save("dynamics/figures/dynamics_$(Lx)x$(Ly)_$(Nt)_$(τ)_examples.png",fig)
