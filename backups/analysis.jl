using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")

dataname = "analysis/data"

Lx = 8
Ly = 8
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 2
T = 0.01

α = 0.9
Nt = 1000
N = 1000

@load "$(dataname)/dara_SA_$(Lx)x$(Ly)_$(T0)_$(T)_$(α)_$(Nt)_$(N).jld2" data
lsT = data["T"]
lsE0 = data["E"]
lsS0 = data["S"]

lsE = zeros(length(lsT))
lsC = zeros(length(lsT))
lsM = zeros(length(lsT))
lsχ = zeros(length(lsT))
for (i,T) in enumerate(lsT)
    lsE[i] = mean(lsE0[i])
    lsC[i] = (mean(lsE0[i] .^ 2) - mean(lsE0[i])^2) / T^2
    lsM[i] = mean(norm.(mean.(lsS0[i])))
    lsχ[i] = (mean(norm.(mean.(lsS0[i])) .^ 2) - mean(norm.(mean.(lsS0[i]))) ^ 2) / T
end

fig = Figure()
figsize = (height = 100,width = 300)
axE = Axis(fig[1,1];xscale = log10,figsize...,
ylabel = L"E")
axM = Axis(fig[2,1];xscale = log10,figsize...,
ylabel = L"M")
axC = Axis(fig[3,1];xscale = log10,figsize...,
ylabel = L"C")
axχ = Axis(fig[4,1];xscale = log10,figsize...,
ylabel = L"\chi")

scatterlines!(axE,lsT,lsE / length(Latt))
scatterlines!(axM,lsT,lsM / length(Latt))
scatterlines!(axC,lsT,lsC / length(Latt))
scatterlines!(axχ,lsT,lsχ / length(Latt))

for ax in [axE,axM,axχ,axC]
    xlims!(ax,extrema(lsT))
end

hidexdecorations!(axE,grid = false,ticks = false)
hidexdecorations!(axM,grid = false,ticks = false)
hidexdecorations!(axC,grid = false,ticks = false)

Label(fig[0,1],text = "$(Lx)x$(Ly) PBC Squa, FM Heisenberg")

resize_to_layout!(fig)
display(fig)

save("analysis/figures/thermal quantity_$(Lx)x$(Ly)_$(T0)_$(T)_$(α)_$(Nt)_$(N).pdf",fig)
save("analysis/figures/thermal quantity_$(Lx)x$(Ly)_$(T0)_$(T)_$(α)_$(Nt)_$(N).png",fig)

