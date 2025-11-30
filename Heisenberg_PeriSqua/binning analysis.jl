
include("../src/src.jl")
include("model.jl")

dataname = "Heisenberg_PeriSqua/data"

Lx = 6
Ly = 6
@load "$(dataname)/Latt_$(Lx)x$(Ly)" Latt

T0 = 2
T = 0.01
params = (J = -1.0,)

α = 0.9
Nt = 1000
N = 1000
W = 5
@load "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT

lsτ = zeros(length(lsT))

for (i,Ttmp) in enumerate(lsT)
    @load "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(round(Ttmp;digits = 8)).jld2" data
    lsW = 2 .^ (0:1:floor(Int64,log(2,N))-5)
    lsK = div.(N, lsW)
    lsτ[i] = BSA(data["E"],lsW,lsK) * W
end
fig = Figure()
ax = Axis(fig[1,1];xscale = log10,
title = "$(Lx)x$(Ly) PBC Squa",
width = 400, height = 250,
xlabel = L"T",
ylabel = L"\tau")

scatterlines!(ax,lsT,lsτ)
lines!(ax,collect(extrema(lsT)), W * ones(2))

xlims!(ax,extrema(lsT))
# ylims!(ax,-0.1)

resize_to_layout!(fig)
display(fig)
save("Heisenberg_PeriSqua/figures/binning analysis_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).pdf",fig)
save("Heisenberg_PeriSqua/figures/binning analysis_$(Lx)x$(Ly)_$(params)_$(T0)_$(T)_$(α)_$(Nt)_$(N).png",fig)
