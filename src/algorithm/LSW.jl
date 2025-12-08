displacement(Latt::AbstractLattice, i::Int64, j::Int64, v::Vector) = (coordinate(Latt,i) - coordinate(Latt,j) + Latt.unitcell.lattice_vecs * v) 

function _LSW_AB(k::Vector,ψ::SimpleState,H::Hamiltonian,vs::Vector;S::Number = 1/2)
    A = zeros(ComplexF64,length(ψ),length(ψ))
    B = zeros(ComplexF64,length(ψ),length(ψ))

    for iname in eachindex(H.name)
        if ndims(H.param[iname]) == 2
            for ((i,j),v) in H.site[iname]
                if i ≠ j
                    A[i,j] += conj(vs[i]') * H.param[iname] * conj(vs[j]) * (S/2) * exp(1im * dot(k, displacement(Latt,i,j,v))) / 2
                    A[j,i] += conj(vs[i]') * H.param[iname] * conj(vs[j]) * (S/2) * exp(-1im * dot(k, displacement(Latt,i,j,v))) / 2
                    B[i,j] += vs[i]' * H.param[iname] * conj(vs[j]) * (S/2) * exp(1im * dot(k, displacement(Latt,i,j,v))) / 2
                    B[j,i] += vs[i]' * H.param[iname] * conj(vs[j]) * (S/2) * exp(-1im * dot(k, displacement(Latt,i,j,v))) / 2 
                else
                    A[i,i] += real(conj(vs[i]') * H.param[iname] * conj(vs[j]) * (S/2) * exp(1im * dot(k, displacement(Latt,i,j,v))))
                    B[i,i] += real(vs[i]' * H.param[iname] * conj(vs[j]) * (S/2) * exp(1im * dot(k, displacement(Latt,i,j,v))))
                end
            end
            for i in 1:length(ψ)
                for j in H.relative[iname][i]
                    if j == i 
                        A[i,i] += - sum(ψ[i]' * H.param[iname] * ψ[j]) * S^2 * 2
                    else
                        A[i,i] += - sum(ψ[i]' * H.param[iname] * ψ[j]) * S^2
                    end
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
    As = zeros(length(lsk))
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
        # @timeit to "build AB" 
        Ak,Bk = _LSW_AB(k,ψ,H,vs;S = S)
        # @timeit to "build AB" 
        Amk,_ = _LSW_AB(-k,ψ,H,vs;S = S)
        # @timeit to "eigen" 
        f = eigen(vcat(hcat(Ak,Bk),-hcat(Bk',conj(Amk))))
        band[:,ik] = f.values[end - length(ψ) + 1:end]
        isweight && (vecs[:,:,ik] = f.vectors[:,end - length(ψ) + 1:end])
        # if mod(ik,showperstep) == 0
        #     show(to;title = "$(ik)/$(length(lsk))")
        #     print("\n")
        # end
        As[ik] = tr(Ak)
    end
    
    @assert mean(abs.(imag.(band))) < 1e-8 "spin not stable"
    band = real.(band)
    ΔE = (sum(band;dims = 1)[:] - As)/2

    # show(to;title = "LSW")
    print("\n")
    
    if isweight
        weight = zeros(length(ψ),length(lsk))
        localaxis = _local_axis.(ψ)
        chiralsp = map(x -> x[1] + 1im * x[2], localaxis)
        chiralsm = map(x -> x[1] - 1im * x[2], localaxis)
        for (ik,k) in enumerate(lsk)
            for i in 1:length(ψ)
                u = vecs[1:length(ψ),i,ik]
                v = vecs[length(ψ)+1:end,i,ik]
                M = zeros(ComplexF64, 3)
                for j in 1:length(ψ)
                    M += exp(1im * dot(k,coordinate(Latt,j))) * (conj(u[j]) * chiralsp[j] + conj(v[j]) * chiralsm[j])
                end
                weight[i,ik] = norm(k) ≈ 0 ? norm(M) ^ 2 : norm(M) ^ 2 - abs2(dot(M[1:2],k)) / norm(k)^2
                # abs2(sum([exp(1im * dot(k,coordinate(Latt,j))) * (u[j] + v[j]) for j in 1:length(ψ)]))
            end
        end
        return band,ΔE,weight
    else
        return band,ΔE
    end
end
