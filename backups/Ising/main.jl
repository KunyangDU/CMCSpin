using JLD2,LatticeUtilities,CairoMakie
include("../lattice/lattice.jl")
include("../src/src.jl")
include("model.jl")

Lx = 6
Ly = 6
Latt = PeriSqua(Lx,Ly)
L = Lx*Ly
nb = neighbor(Latt)
lsβ = 0.1:0.05:0.7
for β in lsβ
times = 50
batch =3000
params = (J=-1.,h=0.)

mag,states = let 
    s = rand(0:2^L-1)
    mag = zeros(batch)
    states = zeros(Int64,batch)
    for i in 1:times
        for j in 1:batch
            ind = rand(1:length(Latt))
            s′ = flip(s,ind)
            accprob(s1,s2) = acceptprob(β,x -> energy(Latt,x,ind,nb;params...),s1,s2)
            isaccept(accprob,s,s′) && (s = s′)
            mag[j] = 2*count1s(s) - L
            states[j] = s
        end
    end
    mag,states
end

magcounts = [sum(map(x -> x == i,mag))/batch for i in -L:2:L]

figsize = (height = 150,width = 150)

finalstate = sum(map(x -> [c - '0' for c in string(x,base = 2,pad = L)], states)) ./ batch
x,y = eachcol(hcat([coordinate(Latt,i) for i in 1:L]...)')
fig = Figure()
axm = Axis(fig[1,1:2];width = figsize.width*2.45,height = figsize.height,
title = "Classical Monte Carlo Simulation\n2D Ising Model, β = $(β)",
xlabel = L"\mathrm{Step}\ n",ylabel = L"\langle m\rangle",
xticks = 0:1000:batch)
lines!(axm,mag ./ L)
ylims!(-1,1)
axfs = Axis(fig[2,1];figsize...,
xticks = (0:L-1,string.(1:L)),yticks = (0:L-1,string.(1:L)),
xlabel = L"x",ylabel = L"y")
heatmap!(axfs,x,y,2finalstate .- 1,colormap = :bwr,colorrange = (-1,1))
axmc = Axis(fig[2,2];figsize...,
xlabel = L"\langle m\rangle",ylabel = L"P(\langle m\rangle)")
barplot!(axmc,(-L:2:L) ./ L,magcounts)
resize_to_layout!(fig)
display(fig)

save("Ising/figures/total_params=$(params)_β=$(β).png",fig)
save("Ising/figures/total_params=$(params)_β=$(β).pdf",fig)
end
