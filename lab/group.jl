include("../src/src.jl")
# include("abstract type.jl")

# using LoopVectorization


function get_heff(ψ::Dict,H::Dict,i::Int64)
    heff = zeros(3)
    for j in eachindex(H["name"])
        if ndims(H["param"][j]) == 1
            i ∈ H["site"][j] && (heff += H["param"][j])
        elseif ndims(H["param"][j]) == 2
            heff += - H["param"][j] * sum(ψ["S"][:,H["relative"][j][i]],dims = 2)[:]
        end
    end
    return heff
end


Lx = 40
Ly = 40

Latt = PeriSqua(Lx,Ly)
# Latt = XCPeriHoneycomb(Lx,Ly)
# Latt = XCPeriTria(Lx,Ly)
groups = group(Latt,[1,],2)
Latt = PeriSqua(Lx,Ly,reverse_map(Tuple(vcat(groups...))))
groups = map(g -> Tuple(Latt.sitemap[g]),groups)
Latt.group = Tuple(groups)

params = (J = 1.0,)

H = Dict(
    "name" => ["J1","h","hsb"],
    "param" => [params.J * diagm([1,1,1]), -[0.0,0.0,0.0], -[0,0,1]],
    "site" => [Tuple(neighbor(Latt)), Tuple(1:length(Latt)), (1,)],
    "relative" => [Tuple(collect(map(x -> x[1] == i ? x[2] : x[1], filter(x -> i in x,Tuple(neighbor(Latt)))))  for i in 1:length(Latt)),nothing,nothing],
    "group" => Latt.group
)

thmalgo = Dict(
    "β" => 1.0,
    "N" => 10,
    "W" => 30,
    "ϵ" => 0.5
)

ψ = Dict(
    "S" => randn(3,length(Latt))
)

normalize!(ψ)


Nthr = get_num_threads_julia()
@time for sites in H["group"]
    L = length(sites)
    Threads.@sync for i in 1:Nthr
        Threads.@spawn begin
            for j in get_chunk_range(L,Nthr,i)
                b = sites[j]
                heff = get_heff(ψ,H,b)
                updata_S_Heatbath(heff,thmalgo["β"])
            end
        end
    end
end

# @time for i in 1:length(Latt)
#     heff = zeros(3)
#     heff += - J * sum(S[:,Jrela[i]],dims = 2)[:]
#     updata_S_Heatbath(heff,β,rs[i]...)
# end

# Jrela[1]

