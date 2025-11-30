
include("../src/src.jl")
include("model.jl")

dataname = "Heisenberg_PeriSqua/data"

Lx = 10
Ly = 10
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 2
T = 0.01
params = (J = -1.0,)

α = 0.9
Nt = 1000
N = 5000
W = 5
@load "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT

lsE = zeros(length(lsT))
lsC = zeros(length(lsT))
lsM = zeros(length(lsT))
lsχ = zeros(length(lsT))
for (i,Ttmp) in enumerate(lsT)
    @load "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(round(Ttmp;digits = 8)).jld2" data
    lsE[i] = mean(data["E"])
    lsC[i] = (mean(data["E"] .^ 2) - mean(data["E"])^2) / Ttmp^2
    lsM[i] = mean(data["M"])
    lsχ[i] = (mean(data["M"] .^ 2) - mean(data["M"])^2) / Ttmp 
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

Label(fig[0,1],text = "$(Lx)x$(Ly) PBC Squa, AFM Heisenberg")

resize_to_layout!(fig)
display(fig)

save("Heisenberg_PeriSqua/figures/thermal quantity_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).pdf",fig)
save("Heisenberg_PeriSqua/figures/thermal quantity_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).png",fig)

