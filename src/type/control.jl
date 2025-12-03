mutable struct ThmAlgo <: AbstractAlgorithm
    W::Int64
    N::Int64
    ϵ::Float64
    function ThmAlgo(W::Int64 = 50, N::Int64 = 1000, ϵ::Float64 = 1.0)
        return new(W,N,ϵ)
    end
end
struct HeatBath <: AbstractAlgorithm end
struct Matropolis <: AbstractAlgorithm
    η::Number
end
struct SDAlgo <: AbstractAlgorithm
    N::Int64
    tol::Float64
    function SDAlgo()
        return new(10000,1e-16)
    end
end

mutable struct UpdateAlgo{S} <: AbstractAlgorithm
    T::Float64
    scheme::S
    function UpdateAlgo(Temp::Float64, scheme::T) where T <: AbstractAlgorithm
        return new{T}(Temp,scheme)
    end
end

struct SAAlgo <: AbstractAlgorithm 
    T0::Float64
    Tf::Float64
    α::Float64
    thermalize::ThmAlgo
    N::Int64
    function SAAlgo(T0::Float64,Tf::Float64,α::Float64, thermalize::ThmAlgo, N::Int64 = 1000)
        return new(T0,Tf,α,thermalize,N)
    end
end