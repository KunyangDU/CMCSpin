
# example:
# - triangular
#     RBASIS3 = [(sqrt(3)/2,1/2,0.),(sqrt(3)/2,-1/2,0.),(0.,0.,1.)]
#     KBASIS2 = kbasis2(RBASIS3)
#     FBZpoint = [(1/6,1/3),(1/3,1/6),(1/6,-1/6),(-1/6,-1/3),(-1/3,-1/6),(-1/6,1/6)]
# - square
#     RBASIS3 = [(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)]
#     KBASIS2 = kbasis2(RBASIS3)
#     FBZpoint = [(1/2,1/2),(1/2,-1/2),(-1/2,-1/2),(-1/2,1/2)]

coordinate(a::Union{Matrix,Vector};basis = KBASIS2) = basism(basis)*a

function kbasis3(basis::Vector)
    basis = collect.(basis)
    V = dot(basis[1],cross(basis[2],basis[3]))
    b1 = cross(basis[1],basis[2])*2*pi/V
    b2 = cross(basis[2],basis[3])*2*pi/V
    b3 = cross(basis[3],basis[1])*2*pi/V
    # return Tuple.([b1,b2,b3])
    return [b1,b2,b3]
end

function kbasis2(basis::Vector)
    kbasis = kbasis3(basis)
    kbasis2 = []
    for kvec in kbasis
        kvec[1] != kvec[2] && push!(kbasis2,kvec[1:2])
    end
    return kbasis2
end

basism(basis::Vector) = hcat(collect.(basis)...)

function isinside(target::Vector,boundary::Matrix;isboundary::Bool = false,tol::Float64 = 1e-8)
    boundaryc = collect.(eachcol(boundary))
    map(x -> push!(x,0),boundaryc)
    targetc = vcat(target,0)
    judge = Vector{Bool}(undef,length(boundaryc))
    judge[end] = true
    std = cross(targetc .- boundaryc[end],boundaryc[1] .- boundaryc[end])[3]
    checkfunc(x,y) = isboundary ? >=(x,y) : >(x,y)
    for i in 1:length(boundaryc)-1
        judge[i] = checkfunc(cross(targetc .- boundaryc[i],boundaryc[i+1] .- boundaryc[i])[3] * std, -tol)
    end
    return sum(judge) == length(boundaryc)
end

isinside(target::Vector,boundary::Vector;basis = KBASIS2,kwargs...) = isinside(target,coordinate(hcat(collect.(boundary)...);basis = basis);kwargs...)
isinside(target::Tuple,boundary::Vector;basis = KBASIS2,kwargs...) = isinside(collect(target),boundary;kwargs...)


v2m(lsv::Vector) = hcat(collect.(lsv)...)

function vrange(lsv::Vector,N::Int64)
    Ls = [norm(lsv[i+1] .- lsv[i]) for i in 1:length(lsv)-1]
    Ns = Int.(round.(N*Ls/sum(Ls);digits=0))
    rnode = vcat(0,[sum(Ls[1:i]) for i in eachindex(Ls)])
    vpath = []
    rpath = []

    for (i,n) in enumerate(Ns)
        ts = range(0,1,n)
        push!(vpath,[lsv[i] + t*(lsv[i+1].-lsv[i]) for t in ts]...)
        push!(rpath,sum(Ls[1:i-1]) .+ ts*Ls[i]...)
    end
    return vpath,rpath,rnode
end
function orientate3(a::Vector,basis = [[0.,-1.],[-sqrt(3)/2,1/2],[sqrt(3)/2,1/2]])
    A = map(x -> dot(a,x),basis)
    return findfirst(x -> x == maximum(A),A),maximum(A)
end

function getxyzbonds(Latt::Union{SimpleLattice{D,S,L}, CompositeLattice{D,S,L}};
    shift = [0,1],
    direction = [[1/2,sqrt(3)/2],[-1/2,sqrt(3)/2],[1,0]],tol=1e-8) where {D,S,L}
    nb = neighbor(Latt)
    _,Ly = S
    return map(direction) do v
        
        filter(x -> let 
            u = coordinate(Latt,x[1]) .- coordinate(Latt,x[2])
            if abs(u[2]) > 1
                u = u .- sign(u[2])*shift*Ly
            end
            abs(u' *v)
        end < tol ,nb)
    end
end