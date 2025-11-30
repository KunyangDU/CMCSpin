function rk4(y::Union{Number, Vector}, dt::Number, f::Function)
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    
    return @. y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
end

