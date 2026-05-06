include("codes/hugget_module.jl")

using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# BACKWARD EGM STEP
#
# Euler: u'(c_t) = β(1+r_{t+1}) E[u'(c_{t+1})]
# Budget t:   c_t  + a' = (1+r_t)   a  + w z  − T_t
# Budget t+1: c_{t+1}   = (1+r_{t+1})a' + w z' − T_{t+1} − σ_{t+1}(a',z')
# ─────────────────────────────────────────────────────────────────────────────
function egm_step_transition(model, σ_next, r_t, r_next, w, T_t, T_next)
    @unpack N_a, N_z, a_vec, z_vec, P_z, β, u_prime, u_prime_inv, a_min, a_max = model

    σ_new     = zeros(N_a, N_z)
    σ_interps = [LinearInterpolation(a_vec, σ_next[:, iz], extrapolation_bc=Line())
                 for iz in 1:N_z]

    for (iz, z) in enumerate(z_vec)
        EMU = zeros(N_a)
        for (ia_next, a_next) in enumerate(a_vec)
            for iz_next in 1:N_z
                a_pp   = σ_interps[iz_next](a_next)
                c_next = (1+r_next)*a_next + w*z_vec[iz_next] - T_next - a_pp
                if c_next > 0
                    EMU[ia_next] += P_z[iz, iz_next] * u_prime(c_next)
                end
            end
        end

        c_endog = zeros(N_a)
        a_endog = zeros(N_a)
        valid   = falses(N_a)

        for (ia_next, a_next) in enumerate(a_vec)
            EMU[ia_next] > 0 || continue
            c_endog[ia_next] = u_prime_inv(β * (1+r_next) * EMU[ia_next])
            a_endog[ia_next] = (c_endog[ia_next] + a_next - w*z + T_t) / (1+r_t)
            valid[ia_next]   = true
        end

        if sum(valid) < 2
            σ_new[:, iz] .= a_min
            continue
        end

        idx  = sortperm(a_endog[valid])
        ae_s = a_endog[valid][idx]
        ap_s = a_vec[valid][idx]
        pol  = LinearInterpolation(ae_s, ap_s, extrapolation_bc=Line())

        for (ia, a) in enumerate(a_vec)
            ap = a < ae_s[1] ? a_min : pol(a)
            ap = clamp(ap, a_min, min((1+r_t)*a + w*z - T_t - 1e-10, a_max))
            σ_new[ia, iz] = ap
        end
    end

    return σ_new
end

# ─────────────────────────────────────────────────────────────────────────────
# BACKWARD PASS: compute σ_path[1..T+1] starting from terminal SS policy σ_T
# ─────────────────────────────────────────────────────────────────────────────
function solve_backward(model, σ_T, r_path, w, T_path, T_len)
    σ_path          = Vector{Matrix{Float64}}(undef, T_len + 1)
    σ_path[T_len+1] = copy(σ_T)

    for t in (T_len-1):-1:0
        σ_path[t+1] = egm_step_transition(
            model, σ_path[t+2],
            r_path[t+1], r_path[t+2],
            w, T_path[t+1], T_path[t+2])
    end

    return σ_path
end

# ─────────────────────────────────────────────────────────────────────────────
# FORWARD PASS: propagate initial distribution λ_0 through policy path
# ─────────────────────────────────────────────────────────────────────────────
function solve_forward(model, λ_0, σ_path, T_len)
    λ_path    = Vector{Matrix{Float64}}(undef, T_len + 1)
    λ_path[1] = copy(λ_0)

    for t in 0:(T_len-1)
        Q           = get_transition_matrix_young(model, σ_path[t+1])
        λ_path[t+2] = reshape(Q' * vec(λ_path[t+1]), model.N_a, model.N_z)
    end

    return λ_path
end

# ─────────────────────────────────────────────────────────────────────────────
# MAIN TRANSITION SOLVER
#
# Given B_path (bond supply path), iterate on r_path until A_t = B_t for all t.
# Update rule: r_t ← damp * A_ss⁻¹(B_t) + (1−damp) * r_t
# where A_ss⁻¹ inverts the pre-computed steady-state asset demand curve.
# ─────────────────────────────────────────────────────────────────────────────
function solve_transition_path(model, ss0, ss1, T_len, B_path;
                                w=1.0, maxiter=500, tol=1e-3, damp=0.7, verbose=true,
                                r_ss_vec=nothing, A_ss_vec=nothing)
    @unpack N_a, N_z, a_vec, z_vec = model

    r_path = collect(range(ss0.r, ss1.r, length=T_len+1))
    T_path = zeros(T_len + 1)   # T = 0 throughout (no lump-sum tax)

    # Build SS A→r inverse interpolant for the update step
    use_inv = !isnothing(r_ss_vec) && !isnothing(A_ss_vec)
    if use_inv
        idx_s  = sortperm(A_ss_vec)
        A_to_r = LinearInterpolation(A_ss_vec[idx_s], r_ss_vec[idx_s],
                                     extrapolation_bc=Line())
    end

    A_path = zeros(T_len + 1)
    σ_path = Vector{Matrix{Float64}}(undef, T_len + 1)
    λ_path = Vector{Matrix{Float64}}(undef, T_len + 1)
    err    = Inf

    for iter in 1:maxiter
        σ_path = solve_backward(model, ss1.σ, r_path, w, T_path, T_len)
        λ_path = solve_forward(model, ss0.λ, σ_path, T_len)

        for t in 0:T_len
            λ_a         = vec(sum(λ_path[t+1], dims=2))
            A_path[t+1] = sum(a_vec .* λ_a)
        end

        excess = A_path .- B_path
        err    = maximum(abs.(excess[2:T_len]))

        if verbose && (iter == 1 || iter % 25 == 0)
            @printf("  iter %3d: max|A−B|=%.5f  r[2]=%.5f  r[T]=%.5f\n",
                    iter, err, r_path[2], r_path[end])
        end

        err < tol && break

        # Update interior r_path using SS inversion A_ss⁻¹(B_t)
        for t in 1:(T_len-1)
            r_new       = use_inv ? A_to_r(B_path[t+1]) :
                                    r_path[t+1] - sign(excess[t+1]) * 2e-4
            r_path[t+1] = damp * r_new + (1-damp) * r_path[t+1]
        end
        r_path[1]       = ss0.r
        r_path[T_len+1] = ss1.r
    end

    err >= tol && verbose &&
        @printf("  Warning: did not converge (max|A−B|=%.5f)\n", err)

    # Aggregate consumption along the path
    C_path = map(0:T_len) do t
        λ_t, σ_t, r_t = λ_path[t+1], σ_path[t+1], r_path[t+1]
        sum(max((1+r_t)*a_vec[ia] + w*z_vec[iz] - σ_t[ia, iz], 0.0) * λ_t[ia, iz]
            for iz in 1:N_z, ia in 1:N_a)
    end

    return (r=r_path, w=fill(w, T_len+1), B=B_path, A=A_path, C=C_path,
            σ=σ_path, λ=λ_path, converged=(err < tol))
end

# ─────────────────────────────────────────────────────────────────────────────
# MAIN EXECUTION
# ─────────────────────────────────────────────────────────────────────────────

w = 1.0
println("=== Huggett Transition Dynamics: s=0 → s_crit ===\n")

# Step 1: asset demand curve (reused for SS inversion during transition)
println("Step 1: Computing asset demand curve...")
model_base = HuggettEGM(T=0.0, s=0.0)
r_natural, r_test_grid, mean_assets_test, negative_idx, r_supply_grid, bond_supply_test =
    asset_test(model_base, -0.06, -0.001; w=w, r_const=0.7, n_neg=46)
demand_neg = mean_assets_test[negative_idx]

# Step 2: r_bar — equilibrium rate when s=0 (no bonds, A(r_bar)=0)
println("Step 2: Finding r_bar...")
r_bar_vec = find_equilibria(HuggettEGM(T=0.0, s=0.0), r_supply_grid, demand_neg,
                             w=w, verbose=false)
r_bar = r_bar_vec[1]
@printf("  r_bar = %.5f\n", r_bar)

# Step 3: maximal deficit and critical rate
println("Step 3: Finding maximal deficit...")
s_critical, r_eq_crit, r_crit, _, _ =
    find_maximal_deficit(r_supply_grid, demand_neg; s_lo=-0.0047, s_hi=-0.0048, tol=1e-9)
@printf("  s_crit = %.8f,  r_crit = %.6f\n\n", s_critical, r_crit)

# Step 4: solve both steady states
println("Step 4: Solving steady states...")
ss0 = solve_welfare(HuggettEGM(T=0.0, s=0.0),       r_bar,  w)
ss1 = solve_welfare(HuggettEGM(T=0.0, s=s_critical), r_crit, w)
@printf("  SS0: r=%.5f,  W=%.6f\n", ss0.r, ss0.W)
@printf("  SS1: r=%.5f,  W=%.6f\n\n", ss1.r, ss1.W)

# Step 5: build B_path — smooth logistic from 0 to B_T
B_T   = s_critical / r_crit
T_len = 300
t_idx = range(0, 1, length=T_len+1)
logit  = @. 1 / (1 + exp(-12 * (t_idx - 0.5)))
logit  = (logit .- logit[1]) ./ (logit[end] - logit[1])
B_path = B_T .* logit

@printf("Step 5: B_path: B_0=%.4f → B_T=%.4f  (T=%d periods)\n\n", B_path[1], B_path[end], T_len)

# Step 6: transition
println("Step 6: Solving transition path...")
trans = solve_transition_path(
    HuggettEGM(T=0.0, s=0.0),
    ss0, ss1, T_len, B_path;
    w=w, damp=0.7, tol=1e-3, verbose=true,
    r_ss_vec=collect(r_supply_grid),
    A_ss_vec=demand_neg
)
println("\nConverged: ", trans.converged)

# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────

t_axis = 0:T_len

# Interest rate path
p_r = plot(t_axis, trans.r,
           xlabel="Period", ylabel=L"r_t",
           title="Interest Rate Transition", legend=:bottomleft, lw=2,
           color=:steelblue, label=L"r_t", guidefont=font(14))
hline!(p_r, [ss0.r], linestyle=:dash, color=:gray,  lw=1, label=L"r_0 = \bar{r}")
hline!(p_r, [ss1.r], linestyle=:dot,  color=:black, lw=1, label=L"r_T = r^*_{crit}")
display(p_r)

# Market clearing: A_t vs B_t
p_mkt = plot(t_axis, [trans.A trans.B],
             xlabel="Period", ylabel="Assets",
             title="Market Clearing along Transition",
             label=["Asset demand " * L"\mathcal{A}_t"
                    "Bond supply "  * L"B_t"],
             lw=2, color=[:steelblue :tomato],
             linestyle=[:solid :dash], guidefont=font(14))
display(p_mkt)

# Aggregate consumption
p_C = plot(t_axis, trans.C,
           xlabel="Period", ylabel=L"C_t",
           title="Aggregate Consumption", legend=false, lw=2,
           color=:seagreen, guidefont=font(14))
display(p_C)

# Asset distribution at selected dates
t_show   = [1, T_len÷4+1, T_len÷2+1, T_len+1]
lbl_show = ["t=0", "t=$(T_len÷4)", "t=$(T_len÷2)", "t=$T_len"]
col_show = [:steelblue, :darkorange, :seagreen, :tomato]

cdf0    = cumsum(vec(sum(trans.λ[1], dims=2)))
idx_pct = something(findfirst(x -> x >= 0.99, cdf0), model_base.N_a)

p_dist = plot(xlabel=L"a", ylabel=L"\mu_a(a)",
              title="Asset Distribution along Transition", guidefont=font(14))
for (t_i, lbl, col) in zip(t_show, lbl_show, col_show)
    λ_a_t = vec(sum(trans.λ[t_i], dims=2))
    plot!(p_dist, model_base.a_vec[1:idx_pct], λ_a_t[1:idx_pct],
          label=lbl, lw=2, color=col)
end
display(p_dist)

# Welfare: aggregate welfare W_t = sum(λ_t .* V_t)
# Compute value functions backward along the transition
println("\nComputing welfare along transition path...")
V_path = Vector{Matrix{Float64}}(undef, T_len + 1)
V_path[T_len+1] = ss1.V
for t in (T_len-1):-1:0
    @unpack N_a, N_z, a_vec, z_vec, P_z, β, u = model_base
    σ_t = trans.σ[t+1]
    r_t = trans.r[t+1]
    V_next = V_path[t+2]
    V_interps = [LinearInterpolation(a_vec, V_next[:, iz], extrapolation_bc=Line())
                 for iz in 1:N_z]
    V_t = zeros(N_a, N_z)
    for iz in 1:N_z, ia in 1:N_a
        c = (1+r_t)*a_vec[ia] + w*z_vec[iz] - σ_t[ia, iz]
        c <= 0 && continue
        ev = sum(P_z[iz, iz2] * V_interps[iz2](σ_t[ia, iz]) for iz2 in 1:N_z)
        V_t[ia, iz] = u(c) + β * ev
    end
    V_path[t+1] = V_t
end

W_path = [sum(trans.λ[t+1] .* V_path[t+1]) for t in 0:T_len]

p_W = plot(t_axis, W_path,
           xlabel="Period", ylabel=L"W_t",
           title="Aggregate Welfare along Transition", legend=false, lw=2,
           color=:purple, guidefont=font(14))
hline!(p_W, [ss0.W], linestyle=:dash, color=:gray, lw=1)
hline!(p_W, [ss1.W], linestyle=:dot,  color=:black, lw=1)
display(p_W)

# CEV: consumption-equivalent welfare change (transition vs staying at ss0)
cev_trans = welfare_cev(W_path[1], ss0.W, model_base.γ)
cev_final = welfare_cev(ss1.W,    ss0.W, model_base.γ)
@printf("\nCEV (start of transition vs ss0) = %.4f%%\n", cev_trans * 100)
@printf("CEV (ss1 vs ss0)                = %.4f%%\n", cev_final * 100)

println("\n=== Done ===")
