using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")


Lx = 6
Ly = 6
Latt = PeriSqua(Lx,Ly)
@load "Heisenberg/data/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
η = 0.02
β = 10
N = 1000000
n = 1000

θ = 0.0 * pi 
ϕ = 0.0 * pi
ĥ = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
@load "Heisenberg/data/total_data_Matropolis_$(Lx)x$(Ly)_$(β)_$(η)_$(N)_$(n)" total_data

Neff = div(N,n)
Ncut = 100
Ni = 10

lsE = zeros(Neff - Ncut)
lsEerr = zeros(Neff - Ncut)
lsM = zeros(Neff - Ncut)
for i in 1:Neff - Ncut
    lsE[i] = total_data[i + Ncut]["E"]
    lsEerr[i] = total_data[i + Ncut]["σE"]
    lsM[i] = sum([dot(ĥ,total_data[i + Ncut]["S"][(j,)])*(-1)^sum(Latt[j][2]) for j in 1:length(Latt)])
end


lsEaccu = cumsum(lsE) ./ collect(1:Neff - Ncut)
lsMaccu = cumsum(lsM) ./ collect(1:Neff - Ncut)

lsEerr = [std(lsE[1:i] / length(Latt)) for i in 1:Neff - Ncut]
lsMerr = [std(lsM[1:i] / length(Latt)) for i in 1:Neff - Ncut]

figsize = (height = 200,width = 400)

fig = Figure()
axE = Axis(fig[1,1];xscale = log10,figsize...,
title = "E = $( round(mean(lsE / length(Latt));digits = 6)), σE =  $(round(std(lsE / length(Latt));digits = 6))")
axM = Axis(fig[2,1];xscale = log10,figsize...,
title = "M = $( round(mean(lsM / length(Latt));digits = 6)), σM =  $(round(std(lsM / length(Latt));digits = 6))")

lines!(axE,Ncut:Neff-1,lsEaccu / length(Latt))
lines!(axM,Ncut:Neff-1,lsMaccu / length(Latt))

errorbars!(axE,Ncut:Ni:Neff - 1,lsEaccu[1:Ni:end] / length(Latt),  lsEerr[1:Ni:end])
errorbars!(axM,Ncut:Ni:Neff - 1,lsMaccu[1:Ni:end] / length(Latt),  lsMerr[1:Ni:end])

for ax in [axE,axM]
    xlims!(ax,Ncut,Neff)
end
# ylims!(axE,-2.,-1.4)
# ylims!(axM,0,1)

resize_to_layout!(fig)
display(fig)
