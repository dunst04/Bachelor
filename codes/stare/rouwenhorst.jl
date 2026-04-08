using QuantEcon, LinearAlgebra
function rouwenhorst1(M, ρ, σ, μ=0)

    p = (1 + ρ)/2

    if M == 1
        P = ones(1,1)

    elseif M == 2
        P = [p 1-p;
            1-p p] 
    else
        P = [p 1-p;
            1-p p]
        for i in 3:M
            P_old = P

            P = p .* [P_old zeros(i-1, 1); zeros(1, i-1) 0] +
            (1-p) .* [zeros(i-1, 1) P_old; 0 zeros(1, i-1)] +
            (1-p) .* [zeros(1, i-1) 0; P_old zeros(i-1, 1)] +
            p .* [0 zeros(1, i-1); zeros(i-1, 1) P_old]
            #normalizing rows
            P[2:end-1, :] .= P[2:end-1, :] ./ 2
        end

        
    end

    ω = σ * sqrt((M - 1)/(1 - ρ^2))

    y = collect(LinRange(μ-ω, μ+ω, M))

    return (; P, y)
end

sol1 = rouwenhorst1(5, 0.9, 0.02)
sol2 = rouwenhorst(5, 0.9, 0.02)

#comparison
maximum(abs.(sol1.P - sol2.p)) 
maximum(abs.(sol1.y - sol2.state_values))