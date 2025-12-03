
function setgroup!(H::Hamiltonian, group::Tuple)
    H.group = group
    for rel in H.relative
        isnothing(rel) && continue
        for g in H.group,i in g 
            @assert !issubset(rel[i], g) "group divided wrong"
        end
    end
    return H
end

function measure(ψ::SimpleState, H::Hamiltonian)
    E = 0.0
    for i in eachindex(H.name)
        if ndims(H.param[i]) == 1
            for s in H.site[i]
                E += H.param[i]' * ψ[s]
            end
        elseif ndims(H.param[i]) == 2
            for (s1,s2) in H.site[i]
                E += ψ[s1]' * H.param[i] * ψ[s2]
            end
        end
    end
    return E
end
