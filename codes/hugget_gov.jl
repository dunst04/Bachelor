### MAIN EXECUTION — Huggett Model with Government

include("hugget_module_gov.jl")

# Create model
model = HuggettEGM()

println("\n=== Huggett Model with Government (EGM) ===\n")
println("Parameters:")
println("  β = $(model.β)")
println("  γ = $(model.γ)")
println("  Borrowing constraint: a ≥ $(model.a_min)")
println("  Asset grid: [$(model.a_min), $(model.a_max)] with $(model.N_a) points")
println("  Income states: $(model.N_z)")
println("\nGovernment:")
println("  Lump-sum tax:       T = $(model.T)")
println("  Government spending: g = $(model.g)")
println("  Permanent deficit:   s = $(model.s)  (= g - T = $(model.g - model.T))")
println("\nProductivity grid: ", round.(model.z_vec, digits=4))

# Solve at initial guess r = -0.01
r_guess = -0.01
w = 1.0

println("\n=== Solving at r = $r_guess ===\n")
@time σ_egm, iter_egm, err_egm = solve_hugget_egm(model, r_guess, w)

# Plot policy function
p_policy = plot_policy_function(model, σ_egm, title_suffix=" (r=$r_guess, gov)")
display(p_policy)

# Compute distributions
λ_egm, _, λ_a_egm, λ_z_egm = stationary_distribution(model, σ_egm)

p_dist_a, p_dist_z = plot_distributions(model, λ_a_egm, λ_z_egm, r_label=" (r=$r_guess)")
display(plot(p_dist_a, p_dist_z, layout=(1,2), size=(1000, 400)))

# Asset demand at r_guess
asset_demand = sum(model.a_vec .* λ_a_egm)
bond_supply   = model.s / r_guess
println("\nAt r = $r_guess:")
println("  Mean assets (demand): $(round(asset_demand, digits=4))")
println("  Bond supply (s/r):    $(round(bond_supply,   digits=4))")
println("  Excess demand:        $(round(asset_demand - bond_supply, digits=4))  (should be ≈ 0 at equilibrium)")

# Euler equation residuals for r_guess
println("\n=== Euler Equation Residuals (r = $r_guess) ===\n")
p_euler_zoom, p_euler_full = plot_euler_residuals(model, σ_egm, r_guess, w,
                                                   title_suffix=" (r=$r_guess)")
display(plot(p_euler_zoom, p_euler_full, layout=(1,2), size=(1400, 500)))

### ASSET SUPPLY vs DEMAND CURVE

println("\n=== Computing Asset Supply vs Demand Curve ===\n")

r_test_grid = range(-0.2, -0.001, length=20)
mean_assets_test = zeros(length(r_test_grid))

println("Computing asset demand for different interest rates...")
for (i, r) in enumerate(r_test_grid)
    print("  r = $(round(r, digits=4))... ")
    σ_temp, _, _ = solve_hugget_egm(model, r, w, verbose=false)
    _, _, λ_a_temp, _ = stationary_distribution(model, σ_temp)
    mean_assets_test[i] = sum(model.a_vec .* λ_a_temp)
    @printf("demand = %.4f  |  supply (s/r) = %.4f\n",
            mean_assets_test[i], -model.s / r)
end

# Bond supply curve: B(r) = -s/r
bond_supply_curve = -model.s ./ r_test_grid

# Find equilibrium: where demand crosses supply
diff = mean_assets_test .- bond_supply_curve
eq_idx = findfirst(i -> diff[i] * diff[i+1] <= 0, 1:length(diff)-1)

if !isnothing(eq_idx)
    r_lo, r_hi = r_test_grid[eq_idx], r_test_grid[eq_idx+1]
    d_lo, d_hi = diff[eq_idx], diff[eq_idx+1]
    r_eq = r_lo - d_lo * (r_hi - r_lo) / (d_hi - d_lo)
    A_eq = mean_assets_test[eq_idx] + (mean_assets_test[eq_idx+1] - mean_assets_test[eq_idx]) *
           (r_eq - r_lo) / (r_hi - r_lo)
    println("\nEquilibrium: r* ≈ $(round(r_eq, digits=5)),  A* ≈ $(round(A_eq, digits=4))")
else
    r_eq, A_eq = NaN, NaN
    println("\nWarning: equilibrium not found in r_test_grid — expand range.")
end

# Plot: demand curve (r on y-axis, assets on x-axis) + supply curve
p_demand = plot(mean_assets_test, collect(r_test_grid),
                ylabel="Interest rate r", xlabel="Mean assets",
                title="Huggett Model with Government\nAsset Market: Supply vs Demand",
                lw=2, marker=:circle, color=:steelblue, legend=:topright,
                label="Asset demand  A(r)")

plot!(p_demand, bond_supply_curve, collect(r_test_grid),
      lw=2, color=:tomato, linestyle=:dash, marker=:diamond,
      label="Bond supply  s/r  (s=$(model.s))")

if !isnan(r_eq)
    scatter!(p_demand, [A_eq], [r_eq],
             markersize=10, color=:black, markershape=:star5,
             label="Equilibrium  r*=$(round(r_eq,digits=4)),  A*=$(round(A_eq,digits=3))")
end

# Reference line r = 1/β - 1
r_natural = 1/model.β - 1
hline!(p_demand, [r_natural], color=:purple, linestyle=:dot, lw=2,
       label="r = 1/β - 1 = $(round(r_natural, digits=4))")

display(p_demand)

#find_equilibrium(model)
