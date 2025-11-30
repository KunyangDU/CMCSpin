using LinearAlgebra
β = 200

# h = randn(3)
h = [1,1,1]
b = β * norm(h)

cθ = log(exp(-b) + rand()*(exp(b) - exp(-b))) / b
sθ = sqrt(1 - cθ^2)
sϕ,cϕ = sincos(2pi * rand())
S′ = [sθ*cϕ,sθ*sϕ,cθ]

ŵ = h / norm(h)
if abs(ŵ[1]) < 0.6 && abs(ŵ[2]) < 0.6
    û = [0.0, -ŵ[3], ŵ[2]]
else
    û = [-ŵ[2], ŵ[1], 0.0]
end
û /= norm(û)
v̂ = cross(ŵ,û)

hcat(û,v̂,ŵ) * S′
