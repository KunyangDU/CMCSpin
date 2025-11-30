function energy(Latt::AbstractLattice,state::Int64;J::Float64 = 1.0, h::Float64 = 0.0)
    EJ = J*sum(
        map(neighbor(Latt)) do pair
            1-2xor(map(x -> bit(state,x),pair)...)
        end
    )
    Eh = -h*(2*count_ones(state)-length(Latt))
    return EJ + Eh
end

function energy(Latt::AbstractLattice,state::Int64,i::Int64,nb::Vector=neighbor(Latt);J::Float64 = 1.0, h::Float64 = 0.0)
    EJ = J*sum(
        map(neighbor(nb,i)) do pair
            1-2xor(map(x -> bit(state,x),pair)...)
        end
    )
    Eh = -h*bit(state,i)
    return EJ + Eh
end