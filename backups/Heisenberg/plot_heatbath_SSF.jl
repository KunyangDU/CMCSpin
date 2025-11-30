using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")



Lx = 6
Ly = 6
Latt = PeriSqua(Lx,Ly)
@load "Heisenberg/data/Latt_$(Lx)x$(Ly)" Latt

hx,hy,hz = 0.0,0.0,0.0
β = 10.0
N = 50000
n = 200

θ = 0.0 * pi 
ϕ = 0.0 * pi
ĥ = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
@load "Heisenberg/data/total_data_Heatbath_$(Lx)x$(Ly)_$(β)_$(N)_$(n)" total_data

Neff = div(N,n)
SM = zeros(length(Latt),length(Latt))
for i in 1:length(Latt), j in i+1:length(Latt)
    SM[i,j] = total_data[end]["SS"][(i,j)]
end
SM += SM'
SM += diagm(ones(length(Latt)))



lskx = 0.999*pi*range(0,1,50)
lsky = 0.999*pi*range(0,1,50)
lsk = filter(x -> isinside(x,FBZpoint;isboundary = true),[[kx,ky] for kx in lskx,ky in lsky][:])
lstk = map(x -> Tuple(x),lsk)
x = map(lsk) do k
    k[1]
end
y = map(lsk) do k
    k[2]
end

SSF = FT2(Latt,SM,lstk)

figsize = (height = 400,width = 400)

fig = Figure()
ax = Axis(fig[1,1];figsize...,autolimitaspect = true)

hm = heatmap!(ax,x,y,SSF)
# boundary!(ax,FBZpoint;color = :black,linewidth = 1.)
# for k in lstk
#     scatter!(ax,k)
# end

Colorbar(fig[1,2],hm)
resize_to_layout!(fig)
display(fig)

# total_data[end]["SS"][(1,18)]
# SM