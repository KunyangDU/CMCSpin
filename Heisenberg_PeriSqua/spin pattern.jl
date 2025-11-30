
include("../src/src.jl")
include("model.jl")

dataname = "Heisenberg_PeriSqua/data/mt"

Lx = 40
Ly = 40
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 2
T = 0.01
params = (J = -1.0,)

α = 0.9
# Tf = T0*(α)^(ceil(Int64,log(0.9,0.01/2)))
Nt = 1000
N = 10000
W = 2
@load "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT
Tf = lsT[end]
@load "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(round(Tf;digits = 8)).jld2" data

lsE = data["E"]
S = data["S"]
# lsE
# S = mean(lsS)

colors = get(colorschemes[:bwr],map(x -> x[3],S),(-1,1))
intensity = 0.8
figsize = (width = 300,height = 300)
fig = Figure()
ax = Axis(fig[1,1];autolimitaspect = true,figsize...,
title = "$(Lx)x$(Ly) PBC Squa")

plotLatt!(ax,Latt)

for i in 1:length(Latt)
    arrowc!(ax,coordinate(Latt,i)...,intensity*S[i][1],intensity*S[i][3];color = colors[i],linewidth = 2)
end

Colorbar(fig[1,2],colorrange = (-1,1),colormap = :bwr,label = L"S_z")

resize_to_layout!(fig)
display(fig)
save("Heisenberg_PeriSqua/figures/spin pattern_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).pdf",fig)
save("Heisenberg_PeriSqua/figures/spin pattern_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).png",fig)

Tf