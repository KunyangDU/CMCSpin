
function SA!(ψ::SimpleState{D,L},H::Hamiltonian,algo::SAAlgo) where {D,L}
    Tmp = algo.T0
    lsE = Float64[]
    lsM = Float64[]
    to = TimerOutput()
    for i in 1:algo.N
        if Tmp < Tf
            show(to;title = "SA $(i-1)/$(ceil(Int64,log(algo.α,algo.Tf/algo.T0)))\nT = $(round(Tmp;digits = 4))")
            return ψ,lsE,lsM
        end
        Tmp = algo.α * Tmp  
        _,tlsE,tlsM,localto,tN = thermalize!(ψ,H,UpdateAlgo(Tmp,HeatBath()),algo.thermalize)
        push!(lsE,tlsE...)
        push!(lsM,tlsM...)
        show(localto;title = "$(i)/$(ceil(Int64,log(algo.α,algo.Tf/algo.T0))) - $(tN)")
        println("\nT=$(Tmp), Tf = $(algo.Tf)")

        merge!(to,localto)
    end
    return ψ,lsE,lsM,to
end
