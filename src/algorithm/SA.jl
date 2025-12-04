
function SA!(ψ::SimpleState{D,L},H::Hamiltonian,algo::SAAlgo;
    showperstep::Int64 = 10, isgc::Bool = false) where {D,L}
    Tmp = algo.T0
    lsE = Float64[]
    lsM = Float64[]
    lsT = Float64[]
    to = TimerOutput()
    for i in 1:algo.N
        push!(lsT,Tmp)
        
        _,tlsE,tlsM,localto,tN = thermalize!(ψ,H,UpdateAlgo(Tmp,HeatBath()),algo.thermalize)
        push!(lsE,tlsE...)
        push!(lsM,tlsM...)
        merge!(to,localto)

        if mod(i,showperstep) == 0
            show(to;title = "$(i)/$(ceil(Int64,log(algo.α,algo.Tf/algo.T0))) - $(tN)")
            println("\nT=$(Tmp), Tf = $(algo.Tf)")
        end

        if Tmp < Tf
            show(to;title = "SA $(i-1)/$(ceil(Int64,log(algo.α,algo.Tf/algo.T0)))\nT = $(round(Tmp;digits = 4))")
            print("\n")
            return ψ,lsT,lsE,lsM
        end
        Tmp = algo.α * Tmp 
        isgc && (@timeit to "manual GC" GC.gc())
    end
    return ψ,lsT,lsE,lsM
end
