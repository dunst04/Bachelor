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

### ASSET DEMAND CURVE
w=1.0
println("\n=== Computing Asset Demand Curve ===\n")

r_natural = 1/model.β - 1
r_test_grid = vcat(range(-0.06, -0.001, length=40),
                   range(0.001, 0.7 * r_natural, length=6))
r_test_grid = collect(r_test_grid)
mean_assets_test = zeros(length(r_test_grid))

println("Computing asset demand for different interest rates...")
for (i, r) in enumerate(r_test_grid)
    print("  r = $(round(r, digits=4))... ")
    σ_temp, _, _ = solve_hugget_egm(model, r, w, verbose=false)
    _, _, λ_a_temp, _ = stationary_distribution(model, σ_temp)
    mean_assets_test[i] = sum(model.a_vec .* λ_a_temp)
    println("mean assets = $(round(mean_assets_test[i], digits=4))")
end

# Bond supply curve B(r) = s/r only exists for r < 0
negative_idx = r_test_grid .< 0
r_supply_grid = r_test_grid[negative_idx]
bond_supply_test = model.s ./ r_supply_grid

# Find all equilibria using the pre-computed demand curve (no re-solving!)
println("\n=== Finding Equilibria ===")
r_equilibria = find_equilibria(model, r_supply_grid, mean_assets_test[negative_idx])

# Plot (x = assets/bonds, y = r)
p_demand = plot(mean_assets_test, collect(r_test_grid),title="s = $(round(model.s, digits=6))",
                ylabel=L"r", xlabel=L"\mathcal{A}",
                lw=2, color=:steelblue, legend=:bottomright, label="asset demand " * L"\mathcal{A}^d(r)", guidefont = font(16), legendfont = font(10))
plot!(p_demand, bond_supply_test, r_supply_grid,
      lw=2, color=:tomato, linestyle=:dash, label="asset supply " * L"\mathcal{A}^s(r) = s/r")

for (k, r_eq) in enumerate(r_equilibria)
    scatter!(p_demand, [model.s / r_eq], [r_eq],
             markersize=4, color=:red, markershape=:circle,
             label="Eq. $k: " * L" r*" * "= $(round(r_eq, digits=4))", guidefont=font(12))
end

hline!(p_demand, [r_natural], color=:purple, linestyle=:dot, lw=2,
       label=L"r = \frac{1}{\beta} - 1" * " = $(round(r_natural, digits=4))")
hline!(p_demand, [0], color=:gray, linestyle=:dot, lw=1, label=false)

display(p_demand)


### MAXIMAL DEFICIT (MINIMAL s) VIA BISECTION

println("\n=== Finding Maximal Deficit (Minimal s) via Bisection ===\n")

# A(r) is independent of s — reuse the pre-computed demand curve.
# find_equilibria uses model.s only for the supply curve B(r) = s/r.
demand_neg = mean_assets_test[negative_idx]

# Bracket: s_lo → 2 equilibria, s_hi → 0 equilibria (given)
s_lo_b, s_hi_b = -0.004, -0.005

n_iter_s = 0
while abs(s_hi_b - s_lo_b) > 1e-6
    s_mid = (s_lo_b + s_hi_b) / 2
    r_eq = find_equilibria(HuggettEGM(s=s_mid), r_supply_grid, demand_neg, verbose=false)
    length(r_eq) == 2 ? (s_lo_b = s_mid) : (s_hi_b = s_mid)
    n_iter_s += 1; n_iter_s >= 200 && break
end

s_critical = (s_lo_b + s_hi_b) / 2
r_eq_crit  = find_equilibria(HuggettEGM(s=s_lo_b), r_supply_grid, demand_neg, verbose=false)
r_crit = mean(r_eq_crit)

println("Converged in $n_iter_s iterations")
@printf("Maximal deficit  -s* =  %.10f\n", -s_critical)
@printf("r*       ≈  %.6f\n\n", mean(r_eq_crit))

# Plot demand and supply for the critical surplus value
bond_supply_critical = s_critical ./ r_supply_grid
p_demand_critical = plot(mean_assets_test, r_test_grid,
                 ylabel=L"r", xlabel=L"\mathcal{A}",
                 title="s = $(round(s_critical, digits=6))",
                 lw=2, color=:steelblue, legend=:bottomright, label="Asset demand " * L"\mathcal{A}^d(r)")
plot!(p_demand_critical, bond_supply_critical, r_supply_grid,
    lw=2, color=:tomato, linestyle=:dash, label="asset supply " * L"\mathcal{A}^s(r) = \frac{s^*}{r}", guidefont=font(16), legendfont=font(10))

for (k, r_eq) in enumerate(r_eq_crit)
    scatter!(p_demand_critical, [s_critical / r_eq], [r_eq],
         markersize=4, color=:red, markershape=:circle,
         label="Critical Eq. $k: " * L"r*" * "= $(round(r_eq, digits=4))")
end

hline!(p_demand_critical, [r_natural], color=:purple, linestyle=:dot, lw=2,
     label=L"r = \frac{1}{\beta} - 1" * " = $(round(r_natural, digits=4))")
hline!(p_demand_critical, [0], color=:gray, linestyle=:dot, lw=1, label=false)

display(p_demand_critical)

solve_and_plot_equilibrium(HuggettEGM(s=s_critical), r_crit, w, verbose=false)





