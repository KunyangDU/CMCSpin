mutable struct Hamiltonian <: AbstractOperator
    name::Vector{String}
    param::Vector{Array}
    site::Vector
    relative::Vector
    group::Union{Nothing,Tuple}
    function Hamiltonian()
        return new(String[], Array[], [], [], nothing)
    end
end


