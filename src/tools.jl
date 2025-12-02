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
