include("../src/src.jl")
using HCubature

Lx = 2
Ly = 2
Latt = PeriSqua(Lx,Ly)

H = Hamiltonian()
addIntr2!(H,Latt,"J1",-diagm([1,1,1]),
Tuple(neighbor_pbc(Latt))
)

ψ = SimpleState(hcat([[0.0,0.0,1.0] * (1) ^ sum(Latt[i][2]) for i in 1:length(Latt)]...))
KB2 = collect(Latt.unitcell.reciprocal_vecs)
# lsk = [KB2 * [i/N,j/N] for i in 0:N-1 for j in 0:N-1]



_minimum_LSWcorrect(ψ,H,Latt;tol = 1e-5)


# I,ϵ = hcubature(x -> _aux_int(ψ,H,Latt,x),(0,0),(1,1),rtol = 1e-6)
# I / Q



