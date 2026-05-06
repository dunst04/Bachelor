### MAIN EXECUTION
# Huggett model without any taxes or transfers
#government issues bonds to finance it's spending
#takes around  7 minutes to run on shitty laptop
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
println(" Lump sum tax: T = $(model.T)")
println(" Government spending: g = $(model.g)")
println(" primary surplus: s = $(model.s)")

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
s_critical, r_eq_crit, r_crit, s_lo_b, s_hi_b = find_maximal_deficit(r_supply_grid, demand_neg, s_lo=-0.0047, s_hi=-0.0048, tol=1e-9)

bond_supply_critical = s_critical ./ r_supply_grid
p_demand_critical = plot_assets(HuggettEGM(s=s_critical), mean_assets_test, r_eq_crit, r_test_grid, r_supply_grid, bond_supply_critical, r_natural, negative_idx)
display(p_demand_critical)

solve_and_plot_equilibrium(HuggettEGM(s=s_critical), r_crit, w, verbose=false)

#welfare analysis for two equilibria s=-0.004
r1, W1, σ1, λ1, λ_a1, λ_z1, V1 = solve_welfare(model, r_equilibria[1], w)
r2, W2, σ2, λ2, λ_a2, λ_z2, V2 = solve_welfare(model, r_equilibria[2], w)

welfare_cev(W2, W1, model.γ)

r0 = find_equilibria(HuggettEGM(T=0.0, g=0.0), r_supply_grid, mean_assets_test[negative_idx])[1]
r3, W3, σ3, λ3, λ_a3, λ_z3, V3 = solve_welfare(HuggettEGM(T=0.0, g=0.0), r0, w)
r4, W4, σ4, λ4, λ_a4, λ_z4, V4 = solve_welfare(HuggettEGM(s=s_critical), r_crit, w)

welfare_cev(W4, W3, model.γ)

p1, p2, p3 =plot_distributions(model, λ_a1, λ_z1)
p4, p5, p6 =plot_distributions(model, λ_a2, λ_z2)
plot(p3, p6)
display(p3)
display(p6)





