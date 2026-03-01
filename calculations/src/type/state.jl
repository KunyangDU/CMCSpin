mutable struct SimpleState{D,L} <: AbstractState
    pattern::Matrix{Float64}
    function SimpleState(A::Matrix{Float64})
        return new{size(A)...}(A)
    end
    function SimpleState(Latt::AbstractLattice, generate::Function = rand)
        A = generate(3,length(Latt))
        return new{size(A)...}(A)
    end
end


function normalize!(ψ::SimpleState{D,L}) where {D,L}
    ψ.pattern ./= norm.(eachcol(ψ.pattern))' .* ones(D)
    return ψ
end

Statistics.mean(ψ::SimpleState) = mean(eachcol(ψ.pattern))
Base.getindex(ψ::SimpleState,i::Int64) = ψ.pattern[:,i]
Base.getindex(ψ::SimpleState,i::Vector{Int64}) = ψ.pattern[:,i]
function Base.setindex!(ψ::SimpleState, a::Vector{Float64}, i::Int64)
    ψ.pattern[:,i] = a
end

Base.length(::SimpleState{D,L}) where {D,L} = L

function Base.iterate(A::SimpleState)
    if isempty(A.pattern)
        return nothing
    else
        return (A.pattern[:,1], 2)
    end
end

function Base.iterate(A::SimpleState, i::Int64)
    if i > length(A)
        return nothing
    else
        return (A.pattern[:,i], i + 1)
    end
end
