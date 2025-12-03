
function addIntr2!(H::Hamiltonian, Latt::AbstractLattice, name::String, param::Matrix, site::Tuple)
    relative = Tuple(collect(map(x -> x[1] == i ? x[2] : x[1], filter(x -> i in x,site)))  for i in 1:length(Latt))
    addIntr!(H,name,param,site,relative)
end
function addIntr1!(H::Hamiltonian, Latt::AbstractLattice, name::String, param::Vector, site::Tuple)
    addIntr!(H,name,param,site,nothing)
end

function addIntr!(H::Hamiltonian,name::String,param::Array,site::Tuple,relative::Union{Nothing,Tuple})
    push!(H.name,name)
    push!(H.param,param)
    push!(H.site,site)
    push!(H.relative,relative)
    return H
end