#initial version of the code was a course material for Quant Econ class by Piotr Zoch
#whom I would like to thank
#source: https://github.com/pzoch/QEcon2025

# Huggett Model with Endogenous Grid Method (EGM)

using Distributions, QuantEcon, Optim, Interpolations, LinearAlgebra, Statistics, ColorSchemes, Plots, Parameters, Printf, LaTeXStrings

@with_kw struct HuggettEGM
    # Preferences
    β = 0.96          # discount factor
    γ = 2.0           # CRRA coefficient
    u = γ == 1 ? (c -> log(c)) : (c -> (c^(1-γ) - 1) / (1-γ))
    u_prime = γ == 1 ? (c -> 1/c) : (c -> c^(-γ))
    u_prime_inv = γ == 1 ? (y -> 1/y) : (y -> y^(-1/γ))
    
    # Income process (Rouwenhorst discretization)
    ρ_z = 0.897                    # persistence
    ν_z = sqrt(0.027)             # volatility
    μ = 0.0                       # mean of log(z)
    N_z = 5                       # number of states
    mc_z = rouwenhorst(N_z, ρ_z, ν_z, μ)
    λ_z = stationary_distributions(mc_z)[1]
    P_z = mc_z.p
    z_vec = exp.(mc_z.state_values) / sum(exp.(mc_z.state_values) .* λ_z)  # normalize mean to 1

    # Asset grid
    ϕ =  0.175              # borrowing constraint
    a_min = -ϕ           # minimum assets
    a_max = 30.0        # maximum assets
    N_a = 500            # grid points

    # Borrowing wedge
    ω_b = 0.5              # extra interest paid by borrowers: r_borrow = r + ω_b

    #government
    T = 0.3 #lump-sum tax
    g = 0.32 #government spending
    s = T - g #government surplus
    
    # Non-uniform grid using polynomial expansion (as in NGM)
    # Concentrates points near a_min where curvature is highest
    θ = 5.0              # curvature parameter (higher = more concentration at lower bound)
    ω = range(0, 1, length=N_a)
    a_vec = a_min .+ (a_max - a_min) .* ω.^θ
end

### EGM ITERATION

function egm_iteration(σ_old, model, r, w)
    """
    One iteration of EGM for given prices (r, w).
    
    Steps:
    1. For each a' on grid, compute expected marginal value using OLD policy
    2. Invert Euler equation to get consumption: c = (u')^{-1}(β(1+r)E[u'(c')])
    3. Back out endogenous asset grid: a = (c + a' - wz)/(1+r)
    4. Interpolate back to exogenous grid
    """
    @unpack N_a, N_z, a_vec, z_vec, P_z, β, u_prime, u_prime_inv, a_min, a_max, T, ω_b = model

    σ_new = zeros(N_a, N_z)

    # Pre-compute interpolations for old policy to speed up the loop
    σ_interps = [LinearInterpolation(a_vec, σ_old[:, iz], extrapolation_bc=Line()) for iz in 1:N_z]

    for (iz, z) in enumerate(z_vec)
        # Step 1: For each a' on grid, compute E[u'(c')] using old policy
        EMU_next = zeros(N_a)  # Expected marginal utility

        for (ia_next, a_next) in enumerate(a_vec)
            # Borrowers (a' < 0) face rate r + ω_b when a' is realized next period
            r_eff_next = a_next < 0 ? r + ω_b : r
            for iz_next in 1:N_z
                z_next = z_vec[iz_next]

                # Use pre-computed interpolation
                a_next_next = σ_interps[iz_next](a_next)

                # Consumption tomorrow
                c_next = (1 + r_eff_next)*a_next + w*z_next - T - a_next_next

                if c_next > 0
                    EMU_next[ia_next] += P_z[iz, iz_next] * u_prime(c_next)
                end
            end
        end

        # Step 2: Invert Euler equation to get consumption on endogenous grid
        # u'(c) = β(1+r_eff(a'))E[u'(c')] => c = (u')^{-1}(β(1+r_eff(a'))E[u'(c')])
        c_endog = zeros(N_a)
        a_endog = zeros(N_a)
        valid = falses(N_a)

        for (ia_next, a_next) in enumerate(a_vec)
            if EMU_next[ia_next] > 0
                r_eff_next = a_next < 0 ? r + ω_b : r

                # Consumption from Euler equation (return on a' depends on sign of a')
                c_endog[ia_next] = u_prime_inv(β * (1 + r_eff_next) * EMU_next[ia_next])

                # Step 3: Back out endogenous asset level today
                # Budget: (1+r_eff(a))*a + w*z = c + a' + T
                # Sign of K determines whether today's a is positive or negative
                K = c_endog[ia_next] + a_next - w*z + T
                r_eff_today = K >= 0 ? r : r + ω_b
                a_endog[ia_next] = K / (1 + r_eff_today)

                valid[ia_next] = true
            end
        end

        # Step 4: Interpolate back to exogenous grid
        if sum(valid) >= 2
            # Sort by endogenous asset grid
            sorted_idx = sortperm(a_endog[valid])
            a_endog_sorted = a_endog[valid][sorted_idx]
            a_next_sorted = a_vec[valid][sorted_idx]

            # Interpolation: given a (exogenous), what is a'?
            policy_interp = LinearInterpolation(a_endog_sorted, a_next_sorted, extrapolation_bc=Line())

            for (ia, a) in enumerate(a_vec)
                # Get policy from interpolation
                a_next_candidate = policy_interp(a)

                # Ensure feasibility: can't save more than total resources
                r_eff_a = a < 0 ? r + ω_b : r
                max_saving = (1 + r_eff_a)*a + w*z - 1e-10  # leave tiny bit for consumption
                a_next_candidate = clamp(a_next_candidate, a_min, min(max_saving, a_max))

                # Handle corner solution at borrowing constraint
                # If current assets 'a' are less than the assets required to choose the
                # lowest a' on the grid (typically a_min), then the constraint binds.
                if a < a_endog_sorted[1]
                    a_next_candidate = a_min
                end

                σ_new[ia, iz] = a_next_candidate
            end
        else
            # Fallback: use old policy if interpolation fails
            #σ_new[:, iz] = σ_old[:, iz]
        end
    end
    
    return σ_new
end

function solve_hugget_egm(model, r, w; maxiter=1000, tol=1e-8, verbose=true)
    """Solve household problem using EGM for given prices."""
    @unpack N_a, N_z, a_vec, z_vec, ω_b = model

    # Initialize policy: save a constant fraction
    σ = zeros(N_a, N_z)
    for (iz, z) in enumerate(z_vec)
        for (ia, a) in enumerate(a_vec)
            r_eff_a = a < 0 ? r + ω_b : r
            income = (1 + r_eff_a)*a + w*z
            σ[ia, iz] = 0.3 * income  # save 30%
            σ[ia, iz] = clamp(σ[ia, iz], a_vec[1], a_vec[end])
        end
    end
    
    err = tol + 1.0
    iter = 1
    
    while err > tol && iter < maxiter
        σ_new = egm_iteration(σ, model, r, w)
        err = maximum(abs.(σ_new - σ))
        σ = σ_new
        iter += 1
    end
    
    if verbose
        println("EGM converged in $iter iterations, error: $err")
    end
    
    return σ, iter, err
end

### STATIONARY DISTRIBUTION

function get_transition_matrix_young(model, σ)
    """
    Compute transition matrix over (a, z) states using Young's method (2010).
    
    Young's method accounts for continuous policy functions by distributing
    probability mass between adjacent grid points through linear interpolation.
    More accurate than discretizing to nearest grid point.
    """
    @unpack N_a, N_z, P_z, a_vec = model
    
    Q = zeros(N_a * N_z, N_a * N_z)
    
    for iz in 1:N_z
        for ia in 1:N_a
            a_next = σ[ia, iz]
            
            # Find grid points that bracket a_next
            if a_next <= a_vec[1]
                # At or below minimum: assign all mass to first point
                ia_next_low = 1
                ia_next_high = 1
                weight_low = 1.0
                weight_high = 0.0
            elseif a_next >= a_vec[end]
                # At or above maximum: assign all mass to last point
                ia_next_low = N_a
                ia_next_high = N_a
                weight_low = 1.0
                weight_high = 0.0
            else
                # Interior: find bracketing points and interpolate
                ia_next_high = findfirst(x -> x >= a_next, a_vec)
                ia_next_low = ia_next_high - 1
                
                # Linear interpolation weights
                a_low = a_vec[ia_next_low]
                a_high = a_vec[ia_next_high]
                weight_high = (a_next - a_low) / (a_high - a_low)
                weight_low = 1.0 - weight_high
            end
            
            # Distribute probability mass according to productivity transitions
            for iz_next in 1:N_z
                # Current state: (ia, iz)
                row = (iz - 1) * N_a + ia
                
                # Next state low: (ia_next_low, iz_next)
                col_low = (iz_next - 1) * N_a + ia_next_low
                Q[row, col_low] += P_z[iz, iz_next] * weight_low
                
                # Next state high: (ia_next_high, iz_next)
                col_high = (iz_next - 1) * N_a + ia_next_high
                Q[row, col_high] += P_z[iz, iz_next] * weight_high
            end
        end
    end
    
    return Q
end

function stationary_distribution(model, σ)
    """
    Compute stationary distribution over (a, z) using Young's method.
    
    Uses matrix power iteration to find the stationary distribution
    of the Markov chain induced by the policy function.
    """
    @unpack N_a, N_z, z_vec = model
    
    Q = get_transition_matrix_young(model, σ)
    
    # Find stationary distribution using power iteration
    # Start with uniform distribution
    λ_vector = ones(N_a * N_z) / (N_a * N_z)
    
    # Iterate until convergence
    for iter in 1:10000
        λ_new = Q' * λ_vector
        
        # Check convergence
        if maximum(abs.(λ_new - λ_vector)) < 1e-10
            λ_vector = λ_new
            break
        end
        
        λ_vector = λ_new
    end
    
    # Normalize (should already be normalized, but ensure numerical stability)
    λ_vector = λ_vector / sum(λ_vector)
    
    # Reshape to (N_a, N_z)
    λ = zeros(N_a, N_z)
    for iz in 1:N_z
        λ[:, iz] = λ_vector[(iz-1)*N_a+1:iz*N_a]
    end
    
    # Compute marginal distributions
    λ_a = sum(λ, dims=2)  # Marginal over assets
    λ_z = sum(λ, dims=1)'  # Marginal over income
    
    return λ, λ_vector, λ_a, λ_z
end

### EQUILIBRIUM

function excess_demand(r, model, w=1.0; verbose=false)
    """
    Compute excess demand for assets at interest rate r.
    Market clears when this equals zero.
    """
    # Solve household problem
    σ, _, _ = solve_hugget_egm(model, r, w, verbose=false)
    
    # Compute stationary distribution
    λ, _, λ_a, _ = stationary_distribution(model, σ)
    
    # Aggregate asset demand
    asset_demand = sum(model.a_vec .* λ_a)
    
    # Bond supply: steady-state government budget  r·B = s  =>  B = s/r
    bond_supply = model.s / r

    #excess demand 
    excess = asset_demand - bond_supply

    if verbose
        @printf("  r = %.5f  A(r) = %.4f  B(r) = s/r = %.4f  excess = %.4f\n",
                r, asset_demand, bond_supply, excess)
    end

    return excess
end

function find_equilibria(model, r_grid, mean_assets; tol=1e-6, w=1.0, verbose=true)
    """
    Find all market-clearing interest rates using a pre-computed asset demand curve.

    Supply curve: B(r) = s/r  (s = T - g; s < 0 → deficit financed by bonds).
    Scans for sign changes of excess demand  A(r) - B(r)  and bisects each interval.
    Returns a vector with 0, 1, or 2 equilibrium interest rates.
    """
    r_grid  = collect(r_grid)
    supply  = model.s ./ r_grid
    excess  = mean_assets .- supply

    r_eq_all = Float64[]
    for i in 1:length(r_grid)-1
        if excess[i] * excess[i+1] < 0
            r_lo, r_hi = r_grid[i], r_grid[i+1]
            e_lo = excess[i]
            for _ in 1:80
                r_mid = (r_lo + r_hi) / 2
                e_mid = excess_demand(r_mid, model, w, verbose=false)
                abs(e_mid) < tol && (r_lo = r_hi = r_mid; break)
                e_lo * e_mid < 0 ? (r_hi = r_mid) : (r_lo = r_mid; e_lo = e_mid)
            end
            push!(r_eq_all, (r_lo + r_hi) / 2)
        end
    end

    if verbose
        if isempty(r_eq_all)
            println("No equilibrium found in scanned range.")
        else
            println("Found $(length(r_eq_all)) equilibrium/equilibria:")
            for (k, r_eq) in enumerate(r_eq_all)
                @printf("  Eq. %d: r* = %.5f   B* = s/r = %.4f\n", k, r_eq, model.s / r_eq)
            end
        end
    end

    return r_eq_all
end

### EULER EQUATION RESIDUALS

function euler_residuals(model, σ, r, w; test_grid=nothing)
    """
    Compute Euler equation residuals to verify solution accuracy.
    
    Euler equation: u'(c) = β(1+r) * E[u'(c')]
    Residual: percentage error in consumption implied by Euler equation
    """
    @unpack a_vec, z_vec, P_z, N_z, β, u_prime, u_prime_inv, T, ω_b = model

    # Create test grid if not provided
    if isnothing(test_grid)
        test_grid = range(model.a_min + 0.01, model.a_max * 0.5, length=500)
    end

    n_test = length(test_grid)
    residuals = zeros(n_test, N_z)

    for (iz, z) in enumerate(z_vec)
        σ_interp = LinearInterpolation(a_vec, σ[:, iz], extrapolation_bc=Line())

        for (ia, a) in enumerate(test_grid)
            a_next = σ_interp(a)
            r_eff_a = a < 0 ? r + ω_b : r
            c = (1 + r_eff_a)*a + w*z - T - a_next

            if c > 0
                r_eff_next = a_next < 0 ? r + ω_b : r
                euler_rhs = 0.0
                for iz_next in 1:N_z
                    z_next = z_vec[iz_next]
                    σ_next_interp = LinearInterpolation(a_vec, σ[:, iz_next], extrapolation_bc=Line())
                    a_next_next = σ_next_interp(a_next)
                    c_next = (1 + r_eff_next)*a_next + w*z_next - T - a_next_next

                    if c_next > 0
                        euler_rhs += P_z[iz, iz_next] * u_prime(c_next)
                    end
                end

                euler_rhs *= β * (1 + r_eff_next)
                c_implied = u_prime_inv(euler_rhs)
                residuals[ia, iz] = (c - c_implied) / c
            else
                residuals[ia, iz] = NaN
            end
        end
    end
    
    return test_grid, residuals
end


### PLOTTING FUNCTIONS

function plot_policy_function(model, σ; a_plot_max=5.0, title_suffix="")
    lines_scheme = get(ColorSchemes.thermal, LinRange(0.2, 0.8, model.N_z))
    
    plot_grid = collect(range(model.a_min, a_plot_max, length=500))

    p = plot(title="Policy Function$title_suffix", xlabel=L"a", ylabel=L"a′", guidefont =font(16))
    for iz in 1:model.N_z
        σ_interp = LinearInterpolation(model.a_vec, σ[:, iz], extrapolation_bc=Line())
        plot!(p, plot_grid, σ_interp.(plot_grid), 
              label="z=$(round(model.z_vec[iz], digits=3))", 
              color=lines_scheme[iz], lw=2)
    end
    plot!(p, plot_grid, plot_grid, label="45°", color=:black, linestyle=:dash)
    
    return p
end

function plot_distributions(model, λ_a, λ_z; r_label="", cdf_pct=0.99, ylim_zoom=0.05)
    """Plot asset and income distributions."""
    p1 = plot(model.a_vec, λ_a, xlabel=L"a", ylabel=L"\mu_a(a)",
            legend=false, lw=2, guidefont=font(16))

    p2 = plot(model.z_vec, λ_z, xlabel=L"z", ylabel=L"\mu_z(z)",
        legend=false, lw=2, linestyle=:dash, marker=:circle, guidefont=font(16))

        # Zoomed asset distribution: x-axis up to cdf_pct percentile, y-axis capped at ylim_zoom
    cdf_a = cumsum(vec(λ_a))
    idx_pct = something(findfirst(x -> x >= cdf_pct, cdf_a), length(cdf_a))
    a_pct = model.a_vec[idx_pct]

    p3 = plot(model.a_vec[1:idx_pct], vec(λ_a)[1:idx_pct],
            xlabel=L"a", ylabel=L"\mu_a(a)",
            title="$(Int(cdf_pct*100))% CDF $r_label",
            legend=false, lw=2, ylims=(0, ylim_zoom), xlims=(model.a_vec[1], a_pct), guidefont=font(16))


    return p1, p2, p3
end

function plot_euler_residuals(model, σ, r, w; test_grid=nothing, title_suffix="")
    """Plot Euler equation residuals with constrained region shading."""
    # Compute residuals
    if isnothing(test_grid)
        test_grid = range(model.a_min, model.a_max, length=2000)
    end
    
    test_a_grid, residuals = euler_residuals(model, σ, r, w, test_grid=test_grid)
    
    # Find y-axis limits
    valid_resids = abs.(residuals[.!isnan.(residuals)])
    ylim_min = minimum(valid_resids[valid_resids .> 0]) / 10
    ylim_max = maximum(valid_resids) * 10
    
    # Select productivity states to plot
    iz_low = 1
    iz_mid = model.N_z ÷ 2 + 1
    iz_high = model.N_z
    iz_indices = [iz_low, iz_mid, iz_high]
    z_labels = ["Low z", "Mid z", "High z"]
    z_colors = [:blue, :red, :green]
    
    # Identify constrained region
    a_constraint_threshold = zeros(model.N_z)
    for iz in 1:model.N_z
        σ_interp = LinearInterpolation(model.a_vec, σ[:, iz], extrapolation_bc=Line())
        for (ia, a) in enumerate(test_a_grid)
            if σ_interp(a) > model.a_min + 1e-6
                a_constraint_threshold[iz] = a
                break
            end
        end
    end
    a_max_constrained = maximum(a_constraint_threshold)
    
    # Panel 1: Zoomed view (constrained region)
    a_zoom_limit = 0.0
    zoom_idx = test_a_grid .<= a_zoom_limit
    
    p_zoom = plot(xlabel="a", ylabel="|Residual|", 
                  title="Euler Residuals: Zoomed$title_suffix",
                  yscale=:log10, ylims=(ylim_min, ylim_max), legend=:topright)
    
    if a_max_constrained > model.a_min
        plot!(p_zoom, [model.a_min, a_max_constrained, a_max_constrained, model.a_min], 
              [ylim_min, ylim_min, ylim_max, ylim_max],
              fillrange=[ylim_min, ylim_min, ylim_min, ylim_min],
              fillalpha=0.2, fillcolor=:gray, linealpha=0, 
              label="Constrained (a'=a_min)")
    end
    
    for (idx, (iz, label, color)) in enumerate(zip(iz_indices, z_labels, z_colors))
        plot!(p_zoom, test_a_grid[zoom_idx], abs.(residuals[zoom_idx, iz]), 
              label="$label (z=$(round(model.z_vec[iz], digits=3)))", 
              linewidth=2, color=color, linestyle=[:solid, :dash, :dot][idx])
    end
    
    # Panel 2: Full view
    p_full = plot(xlabel="a", ylabel="|Residual|", 
                  title="Euler Residuals: Full Grid$title_suffix",
                  yscale=:log10, ylims=(ylim_min, ylim_max), legend=:topright)
    
    if a_max_constrained > model.a_min
        plot!(p_full, [model.a_min, a_max_constrained, a_max_constrained, model.a_min], 
              [ylim_min, ylim_min, ylim_max, ylim_max],
              fillrange=[ylim_min, ylim_min, ylim_min, ylim_min],
              fillalpha=0.2, fillcolor=:gray, linealpha=0, 
              label="Constrained (a'=a_min)")
    end
    
    for (idx, (iz, label, color)) in enumerate(zip(iz_indices, z_labels, z_colors))
        plot!(p_full, test_a_grid, abs.(residuals[:, iz]), 
              label="$label (z=$(round(model.z_vec[iz], digits=3)))", 
              linewidth=2, color=color, linestyle=[:solid, :dash, :dot][idx])
    end
    
    # Print summary statistics
    println("\nMaximum absolute Euler residuals$title_suffix:")
    for (iz, label) in zip(iz_indices, z_labels)
        println("  $label (z=$(round(model.z_vec[iz], digits=3))): ", 
                maximum(abs.(residuals[:, iz])))
    end
    println()
    
    return p_zoom, p_full
end

function solve_and_plot_equilibrium(model, r, w=1.0; verbose=true)
        """
        Solve the household problem for a single interest rate and display
        the policy function, distributions, and Euler residuals in separate figures.
        Returns a named tuple with all computed objects.
        """
        label = "r*=$(round(r, digits=4))"
        verbose && println("\n--- $label ---")

        σ, _, _ = solve_hugget_egm(model, r, w, verbose=verbose)
        λ, _, λ_a, λ_z = stationary_distribution(model, σ)

        p_pol = plot_policy_function(model, σ, title_suffix=" — $label")
        display(p_pol)

        p_dist_a, p_dist_z, p_dist_a_zoom = plot_distributions(model, λ_a, λ_z, r_label=" — $label")
        display(p_dist_a)
        display(p_dist_z)
        display(p_dist_a_zoom)

        p_euler_z, p_euler = plot_euler_residuals(model, σ, r, w, title_suffix=" — $label")
        display(p_euler_z)
        display(p_euler)

        return (r=r, σ=σ, λ_a=λ_a, λ_z=λ_z,
                        p_policy=p_pol, p_dist_a=p_dist_a, p_dist_z=p_dist_z,
                        p_dist_a_zoom=p_dist_a_zoom,
                        p_euler_zoom=p_euler_z, p_euler_full=p_euler)
end
