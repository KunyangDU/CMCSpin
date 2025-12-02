abstract type AbstractLattice end
abstract type SimpleLattice{D,S,L} <: AbstractLattice end

Base.@propagate_inbounds Base.getindex(A::AbstractLattice, i::Int) = reverse(location(Latt,i))
Base.size(A::AbstractLattice) = A.lattice.L
dim(A::AbstractLattice) = dim(A.lattice)
dim(::Lattice{D}) where D = D
Base.length(::SimpleLattice{D,S,L}) where {D,S,L} = L
# Base.@propagate_inbounds Base.setindex!(A::SimpleLattice, value, i::Int) = 1



