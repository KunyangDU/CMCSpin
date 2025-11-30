using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")


Lx = 6
Ly = 6
Latt = PeriSqua(Lx,Ly)
@load "Heisenberg/data/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
η = 0.05
lsβ = 1:20
N = 200000
n = 200

θ = 0.0 * pi 
ϕ = 0.0 * pi
ĥ = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
lsE = zeros(length(lsβ))
lsEerr = zeros(length(lsβ))

for (iβ,β) in enumerate(lsβ)
@show β
@load "Heisenberg/data/total_data_Matropolis_$(Lx)x$(Ly)_$(β)_$(η)_$(N)_$(n)" total_data

Neff = div(N,n)
Ncut = 200
Ni = 50
lsEtmp = zeros(Neff - Ncut)
for i in 1:Neff - Ncut
    lsEtmp[i] += total_data[i + Ncut]["E"]/(Neff - Ncut)
end
lsE[iβ] = mean(lsEtmp)
lsEerr[iβ] = std(lsEtmp)
end
lsE

fig = Figure()

# figsize = (height = 200,width = 400)

# fig = Figure()
ax = Axis(fig[1,1];xscale = log10,figsize...)
scatterlines!(ax,1 ./ lsβ, lsE)


# axM = Axis(fig[2,1];xscale = log10,figsize...,
# title = "M = $( round(mean(lsM / length(Latt));digits = 6)), σM =  $(round(std(lsM / length(Latt));digits = 6))")

# lines!(axE,Ncut:Neff-1,lsEaccu / length(Latt))
# lines!(axM,Ncut:Neff-1,lsMaccu / length(Latt))

# errorbars!(axE,Ncut:Ni:Neff - 1,lsEaccu[1:Ni:end] / length(Latt),  lsEerr[1:Ni:end])
# errorbars!(axM,Ncut:Ni:Neff - 1,lsMaccu[1:Ni:end] / length(Latt),  lsMerr[1:Ni:end])

# for ax in [axE,axM]
#     xlims!(ax,Ncut,Neff)
# end
# # ylims!(axE,-2.,-1.4)
# # ylims!(axM,0,1)

resize_to_layout!(fig)
display(fig)
