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
println("  Lump-sum tax:        T = $(model.T)")
println("  Government spending: g = $(model.g)")
println("  Permanent surplus:   s = T - g = $(model.s)  (s < 0 → deficit)")
println("  Natural rate:        1/β - 1 = $(round(1/model.β - 1, digits=4))")
println("\nProductivity grid: ", round.(model.z_vec, digits=4))

# ── Initial solve at a single r to check policy function ──────────────────────
r_guess = -0.015
w = 1.0

println("\n=== Solving at r = $r_guess ===\n")
@time σ_egm, iter_egm, err_egm = solve_hugget_egm(model, r_guess, w)

# Policy function plot
p_policy = plot_policy_function(model, σ_egm, title_suffix=" (r=$r_guess, gov)")
display(p_policy)

# Distributions
λ_egm, _, λ_a_egm, λ_z_egm = stationary_distribution(model, σ_egm)

p_dist_a, p_dist_z = plot_distributions(model, λ_a_egm, λ_z_egm, r_label=" (r=$r_guess)")
display(plot(p_dist_a, p_dist_z, layout=(1,2), size=(1000, 400)))

# Asset demand vs bond supply at r_guess
asset_demand = sum(model.a_vec .* λ_a_egm)
bond_supply  = model.s / r_guess           # B = s/r  (positive when s<0 and r<0)
println("\nAt r = $r_guess:")
println("  Mean assets (demand):  $(round(asset_demand, digits=4))")
println("  Bond supply (s/r):     $(round(bond_supply,   digits=4))")
println("  Excess demand:         $(round(asset_demand - bond_supply, digits=4))")

# Euler residuals
println("\n=== Euler Equation Residuals (r = $r_guess) ===\n")
p_euler_zoom, p_euler_full = plot_euler_residuals(model, σ_egm, r_guess, w,
                                                   title_suffix=" (r=$r_guess)")
display(plot(p_euler_zoom, p_euler_full, layout=(1,2), size=(1400, 500)))

# ── Supply vs Demand curve + find ALL equilibria ───────────────────────────────

println("\n=== Computing Supply vs Demand Curve and Finding All Equilibria ===\n")

r_nat = 1/model.β - 1

# Scan r < 0: with deficit (s<0), steady-state requires r<0 so that B = s/r > 0
# Both A(r) and B(r)=s/r are upward-sloping for r∈(-0.30,-0.001) → two intersections
r_scan_grid = range(-0.04, -0.014, length=30)
mean_assets  = zeros(length(r_scan_grid))

println("Computing asset demand across interest rates...")
for (i, r) in enumerate(r_scan_grid)
    print("  r = $(round(r, digits=4))... ")
    σ_temp, _, _ = solve_hugget_egm(model, r, w, verbose=false)
    _, _, λ_a_temp, _ = stationary_distribution(model, σ_temp)
    mean_assets[i] = sum(model.a_vec .* λ_a_temp)
    @printf("A(r) = %.4f  |  B(r) = s/r = %.4f\n", mean_assets[i], model.s / r)
end

bond_supply_curve = model.s ./ r_scan_grid    # B(r) = s/r > 0 when s<0 and r<0

# Find ALL equilibria (all sign changes of excess demand)
excess_curve = mean_assets .- bond_supply_curve
r_eq_all = Float64[]
A_eq_all  = Float64[]

for i in 1:length(r_scan_grid)-1
    if excess_curve[i] * excess_curve[i+1] < 0
        # Bisect
        r_lo, r_hi = r_scan_grid[i], r_scan_grid[i+1]
        e_lo = excess_curve[i]
        for _ in 1:80
            r_mid = (r_lo + r_hi) / 2
            e_mid = excess_demand(r_mid, model, w, verbose=false)
            abs(e_mid) < 1e-6 && (r_lo = r_hi = r_mid; break)
            e_lo * e_mid < 0 ? (r_hi = r_mid) : (r_lo = r_mid; e_lo = e_mid)
        end
        r_eq = (r_lo + r_hi) / 2
        push!(r_eq_all, r_eq)
        push!(A_eq_all, model.s / r_eq)
    end
end

println("\n$(length(r_eq_all)) equilibrium/equilibria found:")
for (k, (r_eq, A_eq)) in enumerate(zip(r_eq_all, A_eq_all))
    @printf("  Equilibrium %d:  r* = %.5f   B* = s/r = %.4f\n", k, r_eq, A_eq)
end

# ── Plot: Supply vs Demand with both equilibria ───────────────────────────────

eq_colors  = [:black, :darkgreen, :purple]
eq_shapes  = [:star5, :diamond, :utriangle]

    p_demand = plot(mean_assets, collect(r_scan_grid),
            xlabel="Assets / bonds", ylabel="Interest rate  r",
                title="Huggett with Government — Asset Market\n(surplus s = $(round(model.s,digits=4)), deficit = $(round(-model.s,digits=4)))",
            lw=2, color=:steelblue, label="Asset demand  A(r)",
            xlims=(0, 0.75), legend=:top)

    plot!(p_demand, bond_supply_curve, collect(r_scan_grid),
      lw=2, color=:tomato, linestyle=:dash,
      label="Bond supply  B(r) = s/r  (s = $(model.s))")

for (k, (r_eq, A_eq)) in enumerate(zip(r_eq_all, A_eq_all))
    scatter!(p_demand, [A_eq], [r_eq],
             markersize=11,
             color=eq_colors[mod1(k, length(eq_colors))],
             markershape=eq_shapes[mod1(k, length(eq_shapes))],
             label="Eq. $k:  r* = $(round(r_eq, digits=4)),  B* = $(round(A_eq, digits=3))")
end

    vline!(p_demand, [0.0], color=:gray, linestyle=:dot, lw=1, label="")
    hline!(p_demand, [0.0], color=:black, linestyle=:dot, lw=1, label="r = 0")

display(p_demand)

# ── Solve and show distributions at each equilibrium ─────────────────────────

if !isempty(r_eq_all)
    println("\n=== Solving household problem at each equilibrium ===\n")
    p_policies = []
    p_dists    = []

    for (k, r_eq) in enumerate(r_eq_all)
        println("--- Equilibrium $k: r* = $(round(r_eq, digits=5)) ---")
        σ_eq, iter_eq, err_eq = solve_hugget_egm(model, r_eq, w)
        _, _, λ_a_eq, λ_z_eq  = stationary_distribution(model, σ_eq)

        push!(p_policies, plot_policy_function(model, σ_eq,
              title_suffix=" — Eq. $k (r*=$(round(r_eq,digits=4)))"))
        p_a, _ = plot_distributions(model, λ_a_eq, λ_z_eq,
                                    r_label=" — Eq. $k (r*=$(round(r_eq,digits=4)))")
        push!(p_dists, p_a)
    end

    display(plot(p_policies..., layout=(1, length(r_eq_all)),
                 size=(600*length(r_eq_all), 450)))
    display(plot(p_dists...,    layout=(1, length(r_eq_all)),
                 size=(600*length(r_eq_all), 400)))
end

# ── Nominal equilibria (Taylor rule + ZLB) ────────────────────────────────────

println("\n=== Nominal Equilibria ===")
r_real = isempty(r_eq_all) ? 0.007 : r_eq_all[1]   # use first real equilibrium
eq1_nom, eq2_nom, p_nom = nominal_equilibria(model, r_real)
display(p_nom)

println("\nFiscal Theory interpretation (KNV 2026):")
println("  With s < 0 (deficit) and r* < 0, two real steady-states exist.")
println("  Low-inflation ss (r*_H, B*_H): saddle-path stable (locally determinate).")
println("  High-inflation ss (r*_L, B*_L): locally stable (continuum of equilibria).")
println("  The price level P₀ = B₀/b₀ maps each real ss to a distinct nominal path.")
if eq1_nom[2] >= 0.0
    @printf("  Nom. Eq. 1: π = %+.3f → inflation target,  i = %.3f\n", eq1_nom[1], eq1_nom[2])
else
    @printf("  Nom. Eq. 1: π = %+.3f → target (i = %.3f below ZLB)\n", eq1_nom[1], eq1_nom[2])
end
@printf("  Nom. Eq. 2: π = %+.3f → ZLB trap (i = 0), inflation = |r*|\n", eq2_nom[1])
