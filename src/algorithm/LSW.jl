displacement(Latt::AbstractLattice, i::Int64, j::Int64, v::Vector) = (coordinate(Latt,i) - coordinate(Latt,j) + Latt.unitcell.lattice_vecs * v) / 2

function _LSW_AB(k::Vector,ψ::SimpleState,H::Hamiltonian,vs::Vector;S::Number = 1/2)
    A = zeros(ComplexF64,length(ψ),length(ψ))
    B = zeros(ComplexF64,length(ψ),length(ψ))

    for iname in eachindex(H.name)
        if ndims(H.param[iname]) == 2
            for ((i,j),v) in H.site[iname]
                A[i,j] += conj(vs[i]') * H.param[iname] * conj(vs[j]) * (S/2) * exp(1im * dot(k, displacement(Latt,i,j,v))) / 2
                A[j,i] += conj(vs[i]') * H.param[iname] * conj(vs[j]) * (S/2) * exp(-1im * dot(k, displacement(Latt,i,j,v))) / 2
                B[i,j] += vs[i]' * H.param[iname] * conj(vs[j]) * (S/2) * exp(1im * dot(k, displacement(Latt,i,j,v))) / 2 
                B[j,i] += vs[i]' * H.param[iname] * conj(vs[j]) * (S/2) * exp(-1im * dot(k, displacement(Latt,i,j,v))) / 2 
            end
            for i in 1:length(ψ)
                for j in H.relative[iname][i]
                    A[i,i] += - sum(ψ[i]' * H.param[iname] * ψ[j]) * S^2
                end
            end
        end
    end
    return A,B
end

function LSW(ψ::SimpleState, H::Hamiltonian, lsk::Vector; isweight::Bool = false, S::Number = 1/2,
    showperstep::Int64 = 50)
    vs = map(y -> y[1] + 1im * y[2],map(x -> _local_axis(x),ψ))
    band = zeros(ComplexF64, length(ψ),length(lsk))

    to = TimerOutput()

    if isweight
        vecs = zeros(ComplexF64, 2length(ψ), length(ψ), length(lsk))
    end

    # Nthr = get_num_threads_julia()
    # if Nthr > 1
    #     lsto = [TimerOutput() for _ in 1:Nthr]
    #     counter = Threads.Atomic{Int64}(1)
    #     Threads.@sync for _ in 1:Nthr
    #         # localto = TimerOutput()
    #         Threads.@spawn while true
    #             id = Threads.threadid()
    #             ik = Threads.atomic_add!(counter, 1)
    #             ik > length(lsk) && break
    #             k = lsk[ik]
    #             @timeit lsto[id] "build AB" Ak,Bk = _LSW_AB(k,ψ,H,vs;S = S)
    #             @timeit lsto[id] "build AB" Amk,_ = _LSW_AB(-k,ψ,H,vs;S = S)
    #             @timeit lsto[id] "eigen" f = eigen(vcat(hcat(Ak,Bk),-hcat(Bk',conj(Amk))))
    #             band[:,ik] = f.values[end - length(ψ) + 1:end]
    #             isweight && (vecs[:,:,ik] = f.vectors[:,end - length(ψ) + 1:end])
    #         end
    #     end
    #     for t in lsto
    #         merge!(to, t)
    #     end
    # else
    #     for (ik,k) in enumerate(lsk)
    #         @timeit to "build AB" Ak,Bk = _LSW_AB(k,ψ,H,vs;S = S)
    #         @timeit to "build AB" Amk,_ = _LSW_AB(-k,ψ,H,vs;S = S)
    #         @timeit to "eigen" f = eigen(vcat(hcat(Ak,Bk),-hcat(Bk',conj(Amk))))
    #         band[:,ik] = f.values[end - length(ψ) + 1:end]
    #         isweight && (vecs[:,:,ik] = f.vectors[:,end - length(ψ) + 1:end])
    #     end
    # end

    for (ik,k) in enumerate(lsk)
        @timeit to "build AB" Ak,Bk = _LSW_AB(k,ψ,H,vs;S = S)
        @timeit to "build AB" Amk,_ = _LSW_AB(-k,ψ,H,vs;S = S)
        @timeit to "eigen" f = eigen(vcat(hcat(Ak,Bk),-hcat(Bk',conj(Amk))))
        band[:,ik] = f.values[end - length(ψ) + 1:end]
        isweight && (vecs[:,:,ik] = f.vectors[:,end - length(ψ) + 1:end])
        if mod(ik,showperstep) == 0
            show(to;title = "$(ik)/$(length(lsk))")
            print("\n")
        end
    end

    @assert mean(abs.(imag.(band))) < 1e-8 "spin not stable"
    band = real.(band)

    show(to;title = "LSW")
    print("\n")
    
    if isweight
        weight = zeros(length(ψ),length(lsk))
        for (ik,k) in enumerate(lsk)
            for i in 1:length(ψ)
                u = vecs[1:length(ψ),i,ik]
                v = vecs[length(ψ)+1:end,i,ik]
                weight[i,ik] = abs2(sum([exp(1im * dot(k,coordinate(Latt,j)) / 2) * (u[j] + v[j]) for j in 1:length(ψ)]))
            end
        end
        return band,weight
    else
        return band
    end
end
