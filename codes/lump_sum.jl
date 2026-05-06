#Huggett model with lump sum taxes
include("hugget_module.jl")

model = HuggettEGM(T=0.3, g=0.32)

println("\n=== Huggett Model with Endogenous Grid Method ===\n")
println("Parameters:")
println("  β = $(model.β)")
println("  γ = $(model.γ)")
println("  Borrowing constraint: a ≥ $(model.a_min)")
println("  Asset grid: [$(model.a_min), $(model.a_max)] with $(model.N_a) points")
println("  Income states: $(model.N_z)")
println("\nProductivity grid: ", round.(model.z_vec, digits=4))
println(" Lump sum tax: T = $(model.T)")
println(" Government spending: g = $(model.g)")
println(" primary surplus: s = $(model.s)")

### ASSET DEMAND CURVE
w=1.0
println("\n=== Computing Asset Demand Curve ===\n")

r_natural, r_test_grid, mean_assets_test, negative_idx, r_supply_grid, bond_supply_test =
    asset_test(model, -0.08, -0.005; w=w, r_const=0.7, n_neg=40)

r_equilibria = find_equilibria(model, r_supply_grid, mean_assets_test[negative_idx])

p_demand = plot_assets(model, mean_assets_test, r_equilibria, r_test_grid, r_supply_grid, bond_supply_test, r_natural, negative_idx)
display(p_demand)

### MAXIMAL DEFICIT (MINIMAL s) VIA BISECTION

demand_neg = mean_assets_test[negative_idx]
s_critical, r_eq_crit, r_crit = find_maximal_deficit(r_supply_grid, demand_neg, s_lo = -0.0221850, s_hi = -0.0221851, tol=1e-7)

bond_supply_critical = s_critical ./ r_supply_grid
p_demand_critical = plot_assets(HuggettEGM(s=s_critical), mean_assets_test, r_eq_crit, r_test_grid, r_supply_grid, bond_supply_critical, r_natural, negative_idx)
display(p_demand_critical)

solve_and_plot_equilibrium(HuggettEGM(s=s_critical), r_crit, w, verbose=false)

