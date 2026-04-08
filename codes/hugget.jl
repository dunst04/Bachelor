### MAIN EXECUTION

include("hugget_module.jl")

# Create model
model = HuggettEGM()

println("\n=== Huggett Model with Endogenous Grid Method ===\n")
println("Parameters:")
println("  β = $(model.β)")
println("  γ = $(model.γ)")
println("  Borrowing constraint: a ≥ $(model.a_min)")
println("  Asset grid: [$(model.a_min), $(model.a_max)] with $(model.N_a) points")
println("  Income states: $(model.N_z)")
println("\nProductivity grid: ", round.(model.z_vec, digits=4))

# Solve at initial guess r = 0.01
r_guess = 0.01
w = 1.0

println("\n=== Solving at r = $r_guess ===\n")
@time σ_egm, iter_egm, err_egm = solve_hugget_egm(model, r_guess, w)

# Plot policy function
p_policy = plot_policy_function(model, σ_egm, title_suffix=" (r=$r_guess)")
display(p_policy)

# Compute distributions
λ_egm, _, λ_a_egm, λ_z_egm = stationary_distribution(model, σ_egm)

p_dist_a, p_dist_z = plot_distributions(model, λ_a_egm, λ_z_egm, r_label=" (r=$r_guess)")
display(plot(p_dist_a, p_dist_z, layout=(1,2), size=(1000, 400)))

# Asset demand at r_guess
asset_demand = sum(model.a_vec .* λ_a_egm)
println("\nAt r = $r_guess:")
println("  Mean assets: $(round(asset_demand, digits=4))")
println("  Excess demand: $(round(asset_demand, digits=4)) (should be ≈ 0 at equilibrium)")

# Euler equation residuals for r_guess
println("\n=== Euler Equation Residuals (r = $r_guess) ===\n")
p_euler_zoom, p_euler_full = plot_euler_residuals(model, σ_egm, r_guess, w,
                                                   title_suffix=" (r=$r_guess)")
display(plot(p_euler_zoom, p_euler_full, layout=(1,2), size=(1400, 500)))

### ASSET DEMAND CURVE

println("\n=== Computing Asset Demand Curve ===\n")

r_test_grid = range(-0.2, 0.04, length=15)
mean_assets_test = zeros(length(r_test_grid))

println("Computing asset demand for different interest rates...")
for (i, r) in enumerate(r_test_grid)
    print("  r = $(round(r, digits=4))... ")
    σ_temp, _, _ = solve_hugget_egm(model, r, w, verbose=false)
    _, _, λ_a_temp, _ = stationary_distribution(model, σ_temp)
    mean_assets_test[i] = sum(model.a_vec .* λ_a_temp)
    println("mean assets = $(round(mean_assets_test[i], digits=4))")
end

p_demand = plot(mean_assets_test, r_test_grid, ylabel="Interest rate r", xlabel="Mean assets",
                title="Asset Demand Curve", lw=2, marker=:circle, legend=:topright,
                label="Asset demand")
vline!(p_demand, [0], color=:red, linestyle=:dash, lw=2, label="Market clearing (supply = 0)")

# Add r = 1/β - 1 reference line
r_natural = 1/model.β - 1
hline!(p_demand, [r_natural], color=:purple, linestyle=:dot, lw=2,
       label="r = 1/β - 1 = $(round(r_natural, digits=4))")

display(p_demand)

