
include("../src/src.jl")
include("model.jl")

dataname = "Heisenberg_PeriSqua/data/mt"

Lx = 8
Ly = 8
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 2
T = 0.01
params = (J = 1.0,)

α = 0.9
# Tf = T0*(α)^(ceil(Int64,log(0.9,0.01/2)))
Nt = 1000
N = 10000
W = 2
@load "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT
Tf = lsT[end]
@load "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(N)_$(W)_$(round(Tf;digits = 8)).jld2" data


SS = let SS = zeros(length(Latt),length(Latt))
# lsE = data["E"]
for d in data
    SS += d["SS"] .* length(d["E"]) / N 
end
SS
end


lskx = 0.999*pi*range(0,1,40)
lsky = 0.999*pi*range(0,1,40)
lsk = filter(x -> isinside(x,FBZpoint;isboundary = true),[[kx,ky] for kx in lskx,ky in lsky][:])
lstk = map(x -> Tuple(x),lsk)
x = map(lsk) do k
    k[1]
end
y = map(lsk) do k
    k[2]
end

@time SSF = FT2(Latt,SS,lstk)

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


