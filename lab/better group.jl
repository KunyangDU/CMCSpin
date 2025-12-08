
include("../src/src.jl")

Lx = 60
Ly = 60

Latt = XCPeriTria(Lx,Ly)

# group_optimized(Latt,[1,2,3])
@time groups = group(Latt,[1,2])
# @time groups = group_optimized(Latt,[1,2],3)
@time Latt = XCPeriTria(Lx,Ly,reverse_map(Tuple(vcat(groups...))))
@time groups = map(g -> Tuple(Latt.sitemap[g]),groups)
Latt.group = Tuple(groups)
# groups
fig = Figure()
ax = Axis(fig[1,1],autolimitaspect = true)

plotLatt!(ax,Latt;site = true, tplevel = [1,])

colors = [:red,:blue,:green,:yellow,:purple,:orange,:grey,:black,:gold]
for i in eachindex(groups)
    plotLatt!(ax,Latt;site = true,bond = false,selectedsite = groups[i], sitecolor = repeat([colors[i],],length(groups[i])))
end

resize_to_layout!(fig)
display(fig)

# groups