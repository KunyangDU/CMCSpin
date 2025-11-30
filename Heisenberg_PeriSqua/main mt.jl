# using JLD2,LatticeUtilities,CairoMakie,LinearAlgebra,Statistics,Random
include("../src/src.jl")
include("model.jl")

function calculate!(ψ::Dict,N::Int64,W::Int64,H::Dict,algo::Dict,
    extra!::Function = (x,y) -> nothing)
    to = TimerOutput()
    data = Dict(
        "N" => N,
        "W" => W,
        "E" => Float64[],
        "M" => Float64[],
        "S" => [zeros(3) for _ in 1:length(Latt)]
    )
    @timeit to "calculate" begin
        for j in 1:N
            for k in 1:W 
                @timeit to "update!" update!(ψ,H,algo)
                # merge!(to,local_to;tree_point = ["calculate","update!"])
                data["S"] .+= deepcopy(ψ["S"]) / N / W
                extra!(ψ,data)
            end
            @timeit to "calculate E" push!(data["E"], calculate_E(ψ,H))
            push!(data["M"], norm(sum(ψ["S"])))   
        end
    end

    return ψ,data,to
end


dataname = "Heisenberg_PeriSqua/data/mt"

Lx = 8
Ly = 8
Latt = PeriSqua(Lx,Ly)
@save "$(dataname)/Latt_$(Lx)x$(Ly)" Latt
params = (J = -1.0,)
H = Dict(
    "J1" => (params.J * diagm([1,1,1]),
         Tuple(neighbor(Latt)),
         Dict(i => map(x -> (x[1] == i ? x[2] : x[1]), filter(x -> i in x,Tuple(neighbor(Latt))))  for i in 1:length(Latt))),
    "h" => (-[0.0,0.0,0.0],Tuple(tuple.(1:length(Latt)))),
    "hsb" => (-[0,0,1],((1,),))
)

Obs = Dict(
    "S" => Tuple(tuple.(1:length(Latt))),
    "SS" => (),
)

thmalgo = Dict(
    "N" => 10,
    "W" => 30,
    "ϵ" => 0.5
)

ψ = Dict(
    "S" => [randn(3) |> x -> x/norm(x) for i in 1:length(Latt)]
)

T0 = 2
T = 0.001

α = 0.9
Nt = 1000
N = 100000
W = 2

lsT = Float64[]
let Ttmp = T0
    for i in 1:Nt
        to = TimerOutput()

        if Ttmp < T
            println("T = $(Ttmp) < Tf = $(T) at $(i)th SA")
            break
        end

        Ttmp = α * Ttmp
        algo = Dict(
            "β" => 1 / Ttmp,
            "update" => "Heat bath",
        )

        @timeit to "thermalize" thermalize!(ψ,H,algo,thmalgo)

        Nthr = get_num_threads_julia()
        if Nthr > 1
            Npert = div(N,Nthr)
            data = Dict[]
            tos = TimerOutput[]

            Lock = Threads.ReentrantLock()

            Threads.@sync for _ in 1:Nthr
                Threads.@spawn begin
                    lock(Lock)
                    ψ′ = deepcopy(ψ)
                    unlock(Lock)
                    _,tdata,local_to = calculate!(ψ′,Npert,W,H,algo)
                    tdata["T"] = Ttmp

                    lock(Lock)
                    try
                        push!(data,tdata)
                        push!(tos,local_to)
                    catch
                        rethrow()
                    finally
                        unlock(Lock)
                    end
                end
            end

            if Nthr*Npert ≠ N
                _,fdata,final_to = calculate!(ψ,N- Nthr*Npert,W,H,algo )
                fdata["T"] = Ttmp

                push!(data,fdata)
                push!(tos,final_to)
            end

            for t in tos
                merge!(to,t)
            end

            
            @save "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(N)_$(W)_$(round(Ttmp;digits = 8)).jld2" data
        else
            _,data,local_to = calculate!(ψ,N,W,H,algo)
            merge!(to,local_to)
            @save "$(dataname)/data_$(Lx)x$(Ly)_$(params)_$(N)_$(W)_$(round(Ttmp;digits = 8)).jld2" data
        end

        @timeit to "manual GC" GC.gc()
        push!(lsT,Ttmp)
        show(to;title = "$(i)/$(ceil(Int64,log(α,T/T0))) - T = $(round(Ttmp;digits = 4))")
        print("\n")
    end
end

@save "$(dataname)/lsT_$(T0)_$(T)_$(α)_$(Nt).jld2" lsT


