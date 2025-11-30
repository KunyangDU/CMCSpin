function acceptprob(β::Number,energy::Function,s1::Int64,s2::Int64)
    return exp(-β*(energy(s2) - energy(s1)))
end
function isaccept(accprob::Function,s1::Int64,s2::Int64)
    return rand() < min(accprob(s1,s2),1)
end