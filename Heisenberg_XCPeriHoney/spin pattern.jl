
include("../src/src.jl")
include("model.jl")

dataname = "Heisenberg_XCPeriHoney/data"

Lx = 8
Ly = 8
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 10
T = 0.01
params = (J = 1.0,)

α = 0.9
# Tf = T0*(α)^(ceil(Int64,log(0.9,0.01/2)))
Nt = 1000
N = 1000
W = 5
@load "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT
Tf = lsT[end]
@load "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(round(Tf;digits = 8)).jld2" data

lsE = data["E"]
S = data["S"]
# lsE
# S = mean(lsS)

colors = get(colorschemes[:bwr],map(x -> x[3],S),(-1,1))
intensity = 0.5
figsize = (height = (Ly+1)*80*sqrt(3)/3, width = (Lx + 1)*80)
fig = Figure()
ax = Axis(fig[1,1];autolimitaspect = true,figsize...,
title = "$(Lx)x$(Ly)x2 PBC Honeycomb")

plotLatt!(ax,Latt)

for i in 1:length(Latt)
    arrowc!(ax,coordinate(Latt,i)...,intensity*S[i][1],intensity*S[i][3];color = colors[i],linewidth = 2)
end

Colorbar(fig[1,2],colorrange = (-1,1),colormap = :bwr,label = L"S_z")

resize_to_layout!(fig)
display(fig)
save("Heisenberg_XCPeriHoney/figures/spin pattern_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).pdf",fig)
save("Heisenberg_XCPeriHoney/figures/spin pattern_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).png",fig)

Tf