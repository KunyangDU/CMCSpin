using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")

dataname = "FM Heisenberg/data"

Lx = 10
Ly = 10
Latt = PeriSqua(Lx,Ly)
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
lsT = 0.1:0.1:2.0
N = 100000
n = 100

θ = 0.0 * pi 
ϕ = 0.0 * pi
ĥ = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
lsE = zeros(length(lsT))
lsEerr = zeros(length(lsT))
lsM = zeros(length(lsT))
lsχ = zeros(length(lsT))
lsC = zeros(length(lsT))

for (iT,T) in enumerate(lsT)
@show T
@load "$(dataname)/total_data_Heatbath_$(Lx)x$(Ly)_$(T)_$(N)_$(n)" total_data

Neff = div(N,n)
Ncut = 300
Ni = 20
lsEtmp = zeros(Neff - Ncut)
lsMtmp = zeros(Neff - Ncut)
for i in 1:Neff - Ncut
    lsEtmp[i] = total_data[i + Ncut]["E"]/(Neff - Ncut)
    lsMtmp[i] = norm(sum([total_data[i + Ncut]["S"][(j,)] for j in 1:length(Latt)]))/length(Latt)
end
lsE[iT] = mean(lsEtmp)
lsEerr[iT] = std(cumsum(lsEtmp) ./ (1:Neff - Ncut))
lsM[iT] = mean(lsMtmp)
lsχ[iT] = (mean(lsMtmp .^ 2) - mean(lsMtmp)^2) / T
lsC[iT] = (mean(lsEtmp .^ 2) - mean(lsEtmp)^2) / T^2
end
lsE

fig = Figure()

figsize = (height = 150,width = 300)

# fig = Figure()
axE = Axis(fig[1,1];figsize...,
ylabel = L"E")
axM = Axis(fig[2,1];figsize...,yticks = 0:0.2:1,
ylabel = L"M")
axχ = Axis(fig[3,1];figsize...,
ylabel = L"\chi",xlabel = L"T/J")
# axC = Axis(fig[4,1];figsize...)

scatterlines!(axE,lsT, lsE)
scatterlines!(axM,lsT, lsM)
scatterlines!(axχ,lsT, lsχ)
# scatterlines!(axC,lsT, lsC)

# errorbars!(axE,lsT,lsE,lsEerr)

for ax in [axE,axM,axχ]
    xlims!(ax,extrema(lsT))
end

ylims!(axM,0,1)
ylims!(axχ,0,0.02)

hidexdecorations!(axE,grid = false, ticks = false)
hidexdecorations!(axM,grid = false, ticks = false)

# axM = Axis(fig[2,1];xscale = log10,figsize...,
# title = "M = $( round(mean(lsM / length(Latt));digits = 6)), σM =  $(round(std(lsM / length(Latt));digits = 6))")

# lines!(axE,Ncut:Neff-1,lsEaccu / length(Latt))
# lines!(axM,Ncut:Neff-1,lsMaccu / length(Latt))

# errorbars!(axE,Ncut:Ni:Neff - 1,lsEaccu[1:Ni:end] / length(Latt),  lsEerr[1:Ni:end])
# errorbars!(axM,Ncut:Ni:Neff - 1,lsMaccu[1:Ni:end] / length(Latt),  lsMerr[1:Ni:end])

Label(fig[0,1],text = "$(Lx)x$(Ly) PBC Squa, FM Heisenberg")
# # ylims!(axE,-2.,-1.4)
# # ylims!(axM,0,1)

resize_to_layout!(fig)
display(fig)

save("FM Heisenberg/figures/lsE_$(Lx)x$(Ly)_$(T)_$(N)_$(n).png",fig)
save("FM Heisenberg/figures/lsE_$(Lx)x$(Ly)_$(T)_$(N)_$(n).pdf",fig)
