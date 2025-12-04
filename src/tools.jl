# function FT2(Latt::SimpleLattice, SM::Matrix, k::Tuple)
#     N = length(Latt)
#     @assert SM ≈ SM'
#     # k = [0.0,0.0]
#     SSF = 0.0
#     for i in 1:N, j in i+1:N
#         SSF += SM[i,j] * exp(1im*dot(collect(k),coordinate(Latt,i) - coordinate(Latt,j))) / N
#     end
#     for i in 1:N
#         SSF += SM[i,i] / N
#     end
#     return real(SSF) / N
# end

function FT2(Latt::SimpleLattice, SM::Matrix, k::Tuple)
    N = length(Latt)
    @assert SM ≈ SM'
    Rs = [dot(collect(k) ,coordinate(Latt,i)) for i in 1:N]
    Sqf = exp.(1im * (ones(N) * Rs' .- Rs * ones(N)')) / N
    SSF = sum(Sqf .* SM)
    return real(SSF) / N
end

function FT2(Latt::SimpleLattice, SM::Matrix, tk::Vector)
    return [FT2(Latt,SM,k) for k in tk]
end

function window(x::Number)
    if 1/2 < x ≤ 1
        return -2*(x-1)^3
    elseif 0 ≤ x ≤ 1/2
        return 1 - 6x^2 + 6x^3
    elseif -1/2 ≤ x < 0
        return 1 - 6x^2 - 6x^3
    elseif -1 ≤ x < -1/2
        return 2*(x+1)^3
    else
        return 0
    end
end

get_num_threads_julia() = Threads.nthreads()

function get_chunk_range(N::Int, k::Int, i::Int)
    # i 是当前分块序号 (1 到 k)
    start_idx = div((i - 1) * N, k) + 1
    end_idx   = div(i * N, k)
    return start_idx:end_idx
end

function _local_axis(h::Vector)
    ŵ = h / norm(h)
    if abs(ŵ[1]) < 0.6 && abs(ŵ[2]) < 0.6
        û = [0.0, -ŵ[3], ŵ[2]]
    else
        û = [-ŵ[2], ŵ[1], 0.0]
    end
    û /= norm(û)
    v̂ = cross(ŵ,û)
    return û,v̂,ŵ
end
