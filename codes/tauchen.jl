using LinearAlgebra, SpecialFunctions, QuantEcon

std_norm_cdf(x) = 0.5 * (1 + erf(x / sqrt(2)))

#### APPROXIMATION ####
#### Approximation of AR(1) process using the Tauchen method ####
## x_{t+1} = ρ⋅x_t + ε_t
## ε_t ~ N(μ, σ^2)

function tauchen1(M, ρ, σ, μ=0, m=3)
    
    z = m * σ / sqrt(1- ρ^2)
    y = collect(LinRange(μ-z, μ+z, M))

    D = y[2] - y[1]

    P = Matrix{Float64}(undef, M, M)

    for j in 1:M
        P[j, 1] = std_norm_cdf((y[1] + D/2 - ρ*y[j]) / σ)
        P[j, M] = 1 - std_norm_cdf((y[M] - D/2 - ρ*y[j]) / σ)
        
        for k in 2:M-1
            P[j, k] = std_norm_cdf((y[k] + D/2 - ρ*y[j])/σ) - std_norm_cdf((y[k] - D/2- ρ*y[j])/σ)
        end
    end

    return (; P, y)
end

sol1 = tauchen1(5, 0.9, 0.02)
sol2 = tauchen(5, 0.9, 0.02, 0.0, 3)

#comparison
maximum(abs.(sol1.P - sol2.p))   
maximum(abs.(sol1.y - sol2.state_values))



    

    