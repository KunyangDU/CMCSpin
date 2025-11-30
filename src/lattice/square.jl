mutable struct SquareLattice{D,S,L} <: SimpleLattice{D,S,L}
    unitcell::LatticeUtilities.UnitCell
    lattice::LatticeUtilities.Lattice
    bond::Dict
    # data::Union{Nothing,Vector}
    function SquareLattice(unitcell,lattice::Lattice{D},bond) where D
        S = lattice.L
        return new{D,S,*(S...)}(unitcell,lattice,bond)
    end
end

# function PeriSqua(L,W)
#     square = UnitCell(lattice_vecs = [[1.,0.],[0.,1.]],basis_vecs = [[0.,0.]])
#     lattice = Lattice(L = [W,L], periodic = [true,true])
#     bond = Dict(
#         (true,1) => vcat([[Bond((1,1), [i,0]),Bond((1,1), [0,i])] for i in [-1,1]]...),
#         (false,1) => [Bond((1,1), [1,0]),Bond((1,1), [0,1])],
#         (true,2) => [Bond((1,1), [i,j]) for i in [-1,1],j in [-1,1]][:],
#         (false,2) => [Bond((1,1), [1,i]) for i in [-1,1]],
#     )
#     return SquareLattice(square,lattice,bond)
# end

function PeriSqua(L, W)
    # 正方晶格基矢
    sq = UnitCell(lattice_vecs = [[1., 0.], [0., 1.]], basis_vecs = [[0., 0.]])
    lattice = Lattice(L = [W, L], periodic = [true, true])

    # --- 定义前向偏移量 (Forward Offsets) ---

    # 1近邻 (d=1): 沿轴
    os_1 = [[1, 0], [0, 1]]
    
    # 2近邻 (d=√2): 对角线
    os_2 = [[1, 1], [1, -1]]
    
    # 3近邻 (d=2): 沿轴跳一格 (注意正方晶格没有d=√3的整点)
    os_3 = [[2, 0], [0, 2]]

    # 生成 Bond
    bonds_1_all, bonds_1_fwd = generate_bonds(os_1)
    bonds_2_all, bonds_2_fwd = generate_bonds(os_2)
    bonds_3_all, bonds_3_fwd = generate_bonds(os_3)

    bond = Dict(
        (true, 1) => bonds_1_all, (false, 1) => bonds_1_fwd,
        (true, 2) => bonds_2_all, (false, 2) => bonds_2_fwd,
        (true, 3) => bonds_3_all, (false, 3) => bonds_3_fwd,
    )

    return SquareLattice(sq, lattice, bond)
end
