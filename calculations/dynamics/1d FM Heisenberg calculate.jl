using LatticeUtilities,LinearAlgebra,CairoMakie,LaTeXStrings,JLD2
include("../lattice/lattice.jl")




Lx = 60
Ly = 1
Latt = PeriSqua(Lx,Ly)
params = (J = -1.0,)

nb = filter(x -> x[1] ≠ x[2], neighbor(Latt))
H = Dict(
    "J1" => (params.J * diagm(ones(3)), Tuple(nb)),
    # "h" => (-[hx,hy,hz],Tuple(tuple.(1:length(Latt)))),
    # "hsb" => (-[0.0,0.0,0.0],((1,),))
)


Nt = 1500
τ = 0.02
lst = range(0,Nt*τ,Nt+1)


lsq = range(0, pi, div(Lx,2))
lsω = range(0, 5, div(Lx,1))
# ω = 0.0

dyS = let dyS = zeros(length(lsq),length(lsω))

for i in 1:20
    tmpdyS = zeros(length(lsq),length(lsω))
    ψ = Dict(
        "S" => [[0,0,1] + 0.05*vcat(randn(2),0) |> x -> x/norm(x) for i in 1:length(Latt)]
    )

    lsS = [deepcopy(ψ["S"]),]
    let S = deepcopy(lsS[1])
        for i in 1:Nt 
            S′ = deepcopy(S)
            for j in eachindex(S)

                heff = zeros(3)

                Jij,sites = H["J1"]
                sites = filter(x -> j ∈ x, (sites))
                nbsites = collect(map(x -> j == x[1] ? x[2] : x[1], sites))
                heff += -sum([Jij * S[s] for s in nbsites])

                S′[j] = rk4(S[j], τ, x -> -cross(S[j],heff)) |> x -> x/norm(x)
            end
            push!(lsS, S′)
            S = S′
        end
    end

    for (iq,q) in enumerate(lsq)
        eq = [exp(1im * dot([q,0], coordinate(Latt,i))) for i in 1:length(Latt)] / sqrt(length(Latt))
        Smq0 = sum(eq' .* lsS[1])
        for (iω,ω) in enumerate(lsω)
            lsSq = map( x -> sum(eq .* x), exp.(1im * ω * lst) .* window.(lst/lst[end]) .* lsS) * τ
            tmpdyS[iq,iω] = abs(sum(map(x -> dot(x[1:2],Smq0[1:2]), lsSq)))
        end
    end

    dyS += tmpdyS
end
dyS
end



@save "dynamics/data/dyS_$(Lx)x$(Ly)_$(Nt)_$(τ)_$(params)_examples.jld2" dyS


# fig = Figure()
# ax = Axis(fig[1,1];
# title = "$(Lx)x$(Ly) FM Heisenberg, ϵ=$(round(8/(Nt*τ);digits = 3))",
# width = 300,height = 250,
# xticks = ((0:0.5:1)*pi,[L"0",L"\pi/2",L"\pi"]),
# ylabel = L"\omega",
# xlabel = L"q")

# hm = heatmap!(ax, lsq,lsω, dyS, colorrange = (0,2))
# Colorbar(fig[1,2],hm,label = L"S(q,\omega)")

# xlims!(ax,0,pi)
# ylims!(ax,extrema(lsω))

# resize_to_layout!(fig)

# display(fig)


# save("dynamics/figures/dynamics_FM_$(Lx)x$(Ly)_$(Nt)_$(τ)_examples.pdf",fig)
# save("dynamics/figures/dynamics_FM_$(Lx)x$(Ly)_$(Nt)_$(τ)_examples.png",fig)
