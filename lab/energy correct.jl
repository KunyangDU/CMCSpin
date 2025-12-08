include("../src/src.jl")


Lx = 2
Ly = 2
Latt = PeriSqua(Lx,Ly)

H = Hamiltonian()
addIntr2!(H,Latt,"J1",diagm([1,1,1]),
Tuple(neighbor_pbc(Latt))
)

ψ = SimpleState(hcat([[0.0,0.0,1.0] * (-1) ^ sum(Latt[i][2]) for i in 1:length(Latt)]...))
N = 100
KB2 = collect(Latt.unitcell.reciprocal_vecs)
lsk = [KB2 * [i/N,j/N] for i in 0:N-1 for j in 0:N-1]
_,ΔE = LSW(ψ,H,lsk;isweight = true)
q = _prim_cell(Latt,ψ)[1]
Q = Latt.lattice.N / q

sum(ΔE) / N^2 / Q


