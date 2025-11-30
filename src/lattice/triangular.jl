# 保持你的 Struct 定义不变
mutable struct TriangularLattice{D,S,L} <: SimpleLattice{D,S,L}
    unitcell::LatticeUtilities.UnitCell
    lattice::LatticeUtilities.Lattice
    bond::Dict
    function TriangularLattice(unitcell,lattice::Lattice{D},bond) where D
        S = lattice.L
        return new{D,S,*(S...)}(unitcell,lattice,bond)
    end
end

function XCPeriTria(L, W)
    # 三角晶格基矢
    tria = UnitCell(lattice_vecs = [[1., 0.], [0.5, sqrt(3)/2]], basis_vecs = [[0., 0.]])
    lattice = Lattice(L = [W, L], periodic = [true, true])

    # --- 定义前向偏移量 (Forward Offsets) ---
    
    # 1近邻 (d=1): (1,0), (0,1), (-1,1)
    os_1 = [[1, 0], [0, 1], [-1, 1]]
    
    # 2近邻 (d=√3): (1,1), (1,-2), (2,-1)
    # 这里的几何意义是跨越了两个三角形的高
    os_2 = [[1, 1], [1, -2], [2, -1]]
    
    # 3近邻 (d=2): (2,0), (0,2), (-2,2)
    # 即1近邻的两倍延伸
    os_3 = [[2, 0], [0, 2], [-2, 2]]

    # 生成 Bond
    bonds_1_all, bonds_1_fwd = generate_bonds(os_1)
    bonds_2_all, bonds_2_fwd = generate_bonds(os_2)
    bonds_3_all, bonds_3_fwd = generate_bonds(os_3)

    bond = Dict(
        (true, 1) => bonds_1_all, (false, 1) => bonds_1_fwd,
        (true, 2) => bonds_2_all, (false, 2) => bonds_2_fwd,
        (true, 3) => bonds_3_all, (false, 3) => bonds_3_fwd,
    )

    return TriangularLattice(tria, lattice, bond)
end