### MAIN EXECUTION

include("hugget_module.jl")

# Create model
model = HuggettEGM(T=0.0, g=0.004)

println("\n=== Huggett Model with Endogenous Grid Method ===\n")
println("Parameters:")
println("  β = $(model.β)")
println("  γ = $(model.γ)")
println("  Borrowing constraint: a ≥ $(model.a_min)")
println("  Asset grid: [$(model.a_min), $(model.a_max)] with $(model.N_a) points")
println("  Income states: $(model.N_z)")
println("\nProductivity grid: ", round.(model.z_vec, digits=4))

### ASSET DEMAND CURVE
w=1.0
println("\n=== Computing Asset Demand Curve ===\n")

r_natural, r_test_grid, mean_assets_test, negative_idx, r_supply_grid, bond_supply_test =
    asset_test(model, -0.07, -0.001; w=w, r_const=0.7, n_neg=46)

r_equilibria = find_equilibria(model, r_supply_grid, mean_assets_test[negative_idx])

p_demand = plot_assets(model, mean_assets_test, r_equilibria, r_test_grid, r_supply_grid, bond_supply_test, r_natural, negative_idx)
display(p_demand)

### MAXIMAL DEFICIT (MINIMAL s) VIA BISECTION

demand_neg = mean_assets_test[negative_idx]
s_critical, r_eq_crit, r_crit = find_maximal_deficit(r_supply_grid, demand_neg)

bond_supply_critical = s_critical ./ r_supply_grid
p_demand_critical = plot_assets(HuggettEGM(s=s_critical), mean_assets_test, r_eq_crit, r_test_grid, r_supply_grid, bond_supply_critical, r_natural, negative_idx)
display(p_demand_critical)

solve_and_plot_equilibrium(HuggettEGM(s=s_critical), r_crit, w, verbose=false)





