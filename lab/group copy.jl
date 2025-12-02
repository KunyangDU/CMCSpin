include("../src/src.jl")
# include("model.jl")

Lx = 8
Ly = 8
Latt = XCPeriHoney(Lx,Ly)
# length(Latt)
# Latt.lattice.L
Latt
# @show distance(Latt,1,2)
# sitemap = Tuple(shuffle(1:Lx*Ly))
# asitemap = Tuple([findfirst(x -> x == i,sitemap) for i in 1:Lx*Ly])
# Latt = PeriSqua(Lx,Ly,sitemap)
# @show distance(Latt,1,2)


# nbsites = vcat(neighbor(Latt,1),neighbor(Latt,2))
# initialsites = unique(vcat(collect.(nbsites)...))
# groupsites = Vector[]
# totalsites = Int64[]
# for site in initialsites
#     site in totalsites && continue

#     l = 1
#     sublattsites = [site,]
#     while l ≤ length(sublattsites)
#         nsites = neighborsites(Latt,sublattsites[l];level = 3)
#         for s in nsites
#             if s ∉ sublattsites
#                 push!(sublattsites, s)
#             end
#         end
#         l += 1
#     end
#     sort!(sublattsites)
#     push!(groupsites, sublattsites)
#     push!(totalsites, sublattsites...)

# end

# fig = Figure()
# ax = Axis(fig[1,1])

# plotLatt!(ax,Latt;site = true,sitelabel = true, sitesize = zeros(length(Latt)))
# colors = [:red,:blue,:green,:yellow]
# # for i in eachindex(groupsites)
# #     plotLatt!(ax,Latt;site = true,selectedsite = groupsites[i], sitecolor = repeat([colors[i],],length(groupsites[i])))
# # end

# display(fig)
# groupsites
# groupsites
# neighborsites(Latt,11;level = 2)