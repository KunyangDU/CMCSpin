# using LatticeUtilities

mutable struct HoneycombLattice{D,S,L} <: SimpleLattice{D,S,L}
    unitcell::LatticeUtilities.UnitCell
    lattice::LatticeUtilities.Lattice
    bond::Dict
    function HoneycombLattice(unitcell,lattice::Lattice{D},bond) where D
        S = lattice.L
        return new{D,S,*(S...)}(unitcell,lattice,bond)
    end
end

"""
XCPeriHoney(L, W)

构建双周期（Torus）蜂窝晶格。

参数:
- W: 晶格宽度方向的**元胞数** (沿 a1 方向)
- L: 晶格长度方向的**元胞数** (沿 a2 方向)

总格点数 = 2 * W * L
"""
function XCPeriHoney(L, W)
    # 1. 定义元胞 (Unit Cell)
    # 这里的 lattice_vecs 定义了元胞的形状（三角晶格基底）
    # basis_vecs 定义了元胞内的 2 个格点：A (索引1) 和 B (索引2)
    # 坐标系下: A在原点, B在重心位置 (1/3, 1/3)
    honey = UnitCell(
        lattice_vecs = [[1., 0.], [0.5, sqrt(3)/2]],
        basis_vecs   = [[0., 0.], [0.5, sqrt(3)/6]] 
    )
    
    # 2. 定义晶格整体尺寸 (Lattice Dimensions)
    # L=[W, L] 表示在第一个方向上有 W 个元胞，第二个方向上有 L 个元胞
    # periodic=[true, true] 实现双向周期性边界条件 (Torus)
    lattice = Lattice(L = [W, L], periodic = [true, true])

    # 3. 定义连接 (Bonds)
    # 我们使用 (n1, n2) 整数坐标表示元胞偏移量
    
    # --- 辅助: 三角晶格方向 (用于同子格 A-A, B-B 连接) ---
    tri_fwd = [[1, 0], [0, 1], [-1, 1]]
    tri_all = vcat(tri_fwd, [-1 .* os for os in tri_fwd])

    # --- 1近邻 (Nearest Neighbor, d = 1/√3) ---
    # 连接类型: A(1) -> B(2)
    # 偏移量: [0,0] (元胞内), [-1,0] (左邻元胞), [0,-1] (下邻元胞)
    # 物理距离平方: 1/3
    nn_offsets = [[0, 0], [-1, 0], [0, -1]]
    bonds_1 = [Bond((1, 2), os) for os in nn_offsets]

    # --- 2近邻 (Next Nearest Neighbor, d = 1) ---
    # 连接类型: A(1)->A(1) 和 B(2)->B(2)
    # 实际上就是底层的三角晶格连接
    # 物理距离平方: 1
    bonds_2_fwd = vcat([Bond((1, 1), os) for os in tri_fwd],
                       [Bond((2, 2), os) for os in tri_fwd])
                       
    bonds_2_all = vcat([Bond((1, 1), os) for os in tri_all],
                       [Bond((2, 2), os) for os in tri_all])

    # --- 3近邻 (Third Nearest Neighbor, d = 2/√3) ---
    # 连接类型: A(1) -> B(2)
    # 这些是跨越六边形的对角连接
    # 偏移量: [1, -1], [-1, 1], [-1, -1]
    # 物理距离平方: 4/3
    nnn_offsets = [[1, -1], [-1, 1], [-1, -1]]
    bonds_3 = [Bond((1, 2), os) for os in nnn_offsets]

    bond = Dict(
        # --- Order 1 (1NN) ---
        (true, 1)  => bonds_1,     # 所有连接 (用于几何分析)
        (false, 1) => bonds_1,     # 唯一前向 (A->B 是唯一的，无需剔除反向)

        # --- Order 2 (2NN) ---
        (true, 2)  => bonds_2_all, 
        (false, 2) => bonds_2_fwd, # 同子格连接，只保留前向避免 Hamiltonian 重复计数

        # --- Order 3 (3NN) ---
        (true, 3)  => bonds_3,
        (false, 3) => bonds_3      # 同样是 A->B，唯一前向
    )
    
    return HoneycombLattice(honey, lattice, bond)
end