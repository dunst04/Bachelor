using Plots, LaTeXStrings

# 1. Calibration
β = 0.99        # Discount factor 
κ = 0.1275      # Slope of the Phillips curve

# 2. Define the boundary: ϕ_π = 1 - ((1 - β) / κ) * ϕ_x
# Rearranging the inequality from the user image
taylor_boundary(ϕ_x) = 1 - ((1 - β) / κ) * ϕ_x

# 3. Plotting
ϕ_x_vals = range(0, 4, length=200)
ϕ_π_boundary = taylor_boundary.(ϕ_x_vals)

plot(ϕ_x_vals, ϕ_π_boundary, 
    label="Determinacy Boundary", 
    lw=3, color=:black, 
    xlabel=L"Response to Output Gap $\phi_x$", 
    ylabel=L"Response to Inflation $\phi_\pi$",
    title="Determinacy vs. Indeterminacy Regions",
    xlims=(0, 4), ylims=(0, 2),
    fillrange=2, fillalpha=0.2, fillcolor=:green, 
    legend=false)

# Adding labels for regions
annotate!([(2, 1.5, text("Determinacy\n(Unique Solution)", :green, :center, 10)),
           (0.5, 0.5, text("Indeterminacy", :red, :center, 10))])

# Shading the indeterminacy region below the line
plot!(ϕ_x_vals, ϕ_π_boundary, 
    fillrange=0, fillalpha=0.2, fillcolor=:red, label="")