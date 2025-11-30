function generate_bonds(fwd_offsets)
    # 生成反向偏移量，用于构建 "所有方向" 的连接
    all_offsets = vcat(fwd_offsets, [-1 .* os for os in fwd_offsets])
    
    fwd_bonds = [Bond((1,1), os) for os in fwd_offsets]
    all_bonds = [Bond((1,1), os) for os in all_offsets]
    return all_bonds, fwd_bonds
end

function Base.length(lat::SimpleLattice)
    N_cells = prod(lat.lattice.L)             # W * L
    N_basis = length(lat.unitcell.basis_vecs) # 基底数量 (Tri=1, Honey=2)
    return N_cells * N_basis
end

