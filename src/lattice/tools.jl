function generate_bonds(fwd_offsets, basis_idxs=(1,1))
    # basis_idxs: (from, to)
    # 如果 from != to (如 Honeycomb A->B)，则无须生成反向，全向=前向
    if basis_idxs[1] != basis_idxs[2]
        bonds = [Bond(basis_idxs, os) for os in fwd_offsets]
        return bonds, bonds
    else
        # 同一基底，需要生成反向以获得几何上的“所有邻居”
        all_offsets = vcat(fwd_offsets, [-1 .* os for os in fwd_offsets])
        fwd_bonds = [Bond(basis_idxs, os) for os in fwd_offsets]
        all_bonds = [Bond(basis_idxs, os) for os in all_offsets]
        return all_bonds, fwd_bonds
    end
end
# function Base.length(lat::SimpleLattice)
#     N_cells = prod(lat.lattice.L)             # W * L
#     N_basis = length(lat.unitcell.basis_vecs) # 基底数量 (Tri=1, Honey=2)
#     return N_cells * N_basis
# end

reverse_map(sitemap::Tuple) = Tuple([findfirst(x -> x == i,sitemap) for i in Tuple(1:length(sitemap))])

function group(Latt::AbstractLattice, level0::Vector, level::Int64 = maximum(level0) + 1)
    nbsites = vcat([neighbor(Latt,1;level = i) for i in level0]...)
    initialsites = unique(vcat(collect.(nbsites)...))
    groupsites = Vector[]
    totalsites = Int64[]
    for site in initialsites
        site in totalsites && continue
        l = 1
        sublattsites = [site,]
        while l ≤ length(sublattsites)
            nsites = neighborsites(Latt,sublattsites[l];level = level)
            push!(sublattsites, filter(x -> x ∉ sublattsites, nsites)...)
            l += 1
        end
        sort!(sublattsites)
        push!(groupsites, sublattsites)
        push!(totalsites, sublattsites...)
    end
    @assert isequal(1:length(Latt),sort(totalsites))
    @assert length(groupsites) ≠ 1 "group failed, change the lattice size"
    return groupsites
end

