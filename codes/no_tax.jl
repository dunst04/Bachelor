### MAIN EXECUTION
# Huggett model without any taxes or transfers
#government issues bonds to finance it's spending
#takes around  7 minutes to run on shitty laptop
include("hugget_module.jl")

function plot_both_distributions(model, λ_a1, λ_a2; cdf_pct=0.99,
                                  labels=["Eq. 1", "Eq. 2"], ylim_zoom=0.05)
    cdf1 = cumsum(vec(λ_a1))
    cdf2 = cumsum(vec(λ_a2))
    idx1 = something(findfirst(x -> x >= cdf_pct, cdf1), length(cdf1))
    idx2 = something(findfirst(x -> x >= cdf_pct, cdf2), length(cdf2))
    idx_max = max(idx1, idx2)
    a_max_plot = model.a_vec[idx_max]

    p = plot(model.a_vec[1:idx_max], vec(λ_a1)[1:idx_max],
             xlabel=L"a", ylabel=L"\mu_a(a)",
             title="Asset Distribution ($(Int(cdf_pct*100))% CDF)",
             label=labels[1], lw=2,
             ylims=(0, ylim_zoom), xlims=(model.a_vec[1], a_max_plot),
             guidefont=font(16))
    plot!(p, model.a_vec[1:idx_max], vec(λ_a2)[1:idx_max],
          label=labels[2], lw=2, linestyle=:dash)
    return p
end

function plot_both_policies(model, σ1, σ2; a_plot_max=7.0, labels=["Eq. 1", "Eq. 2"])
    iz_indices = [1, model.N_z ÷ 2 + 1, model.N_z]
    z_labels   = ["Low y", "Mid y", "High y"]
    colors     = [:steelblue, :darkorange, :seagreen]
    plot_grid  = collect(range(model.a_min, a_plot_max, length=500))

    p = plot(title="Policy Functions", xlabel=L"a", ylabel=L"a′", guidefont=font(16))
    for (k, iz) in enumerate(iz_indices)
        σ1_interp = LinearInterpolation(model.a_vec, σ1[:, iz], extrapolation_bc=Line())
        σ2_interp = LinearInterpolation(model.a_vec, σ2[:, iz], extrapolation_bc=Line())
        plot!(p, plot_grid, σ1_interp.(plot_grid), label="$(labels[1]) — $(z_labels[k])", color=colors[k], lw=2)
        plot!(p, plot_grid, σ2_interp.(plot_grid), label="$(labels[2]) — $(z_labels[k])", color=colors[k], lw=2, linestyle=:dash)
    end
    plot!(p, plot_grid, plot_grid, label="45°", color=:black, linestyle=:dot)
    return p
end

# Create model
model = HuggettEGM(T=0.0, s=-0.004)

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
    asset_test(model, -0.06, -0.001; w=w, r_const=0.7, n_neg=46)

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

p_dist_both = plot_both_distributions(model, λ_a1, λ_a2)

p_pol_both = plot_both_policies(model, σ1, σ2)


r_bar = find_equilibria(HuggettEGM(T=0.0, s=0.0), r_supply_grid, mean_assets_test[negative_idx], w=w)
r3, W3, σ3, λ3, λ_a3, λ_z3, V3 = solve_welfare(HuggettEGM(T=0.0, s=s_critical), r_crit, w)
r4, W4, σ4, λ4, λ_a4, λ_z4, V4 = solve_welfare(HuggettEGM(T=0.0, s=0.0), r_bar[1], w)

p_dist_crit = plot_both_distributions(model, λ_a3, λ_a4, labels=["s = $(round(s_critical, digits=5))", "s = 0.0"])

p_pol_crit = plot_both_policies(model, σ3, σ4, labels=["s = $(round(s_critical, digits=5))", "s = 0.0"])

welfare_cev(W3, W4, model.γ)