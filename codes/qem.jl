using LinearAlgebra, Plots

#setting parameters
α=0.33
β=0.99
θ=2
δ=0.025
ϕ=1
ρ=0.95
σ=0.01
Λ=(1-α)/(ϕ+α)
b=1/β-1+δ

n=2
m=1

A1 = [1 0 0; 0 1 0; -β*b*(1+ϕ)/(ϕ+α) β*b*ϕ*Λ θ*(1+β*b*Λ)]
A0 = [ρ 0 0; b/α*(1+Λ) 1/β+b*Λ δ-b/α*(1+θ*Λ); 0 0 θ]
γ = [1; 0; 0]
A = A1 \ A0
R=inv(A1)*γ

decomp = eigen(A)
λ = decomp.values
C = decomp.vectors

p = sortperm(abs.(λ))
C_sort = C[:, p]

Cinv = inv(C_sort)
m_star = sum(abs.(λ) .> 1)
n_star = n + m - m_star
# Checking conditions of Theorem 1
if m_star != m
	println("The Blanchard-Kahn condition is not satisfied.")
end
if det(Cinv[n_star+1:end, n_star+1:end]) == 0
	println("The rank condition is not satisfied.")
end
# Solution
Γ = -inv(Cinv[n_star+1:end, n_star+1:end]) * Cinv[n_star+1:end, 1:n_star]
Φ2= A[1:n_star, 1:n_star] + A[1:n_star, n_star+1:end] * Γ
Ψ2 = R[1:n_star]

Φ = [Φ2 zeros(2,1); Γ*Φ2 0]
Ψ = [Ψ2; Γ*Ψ2]

T=100
a=zeros(T)
k=zeros(T)
c=zeros(T)
y=zeros(T)
l=zeros(T)


W = zeros(3)

for t in 1:T
    if t == 1
       global W = Φ * W + Ψ .* σ
    else 
        global W = Φ * W
    end
    a[t] = W[1]
    k[t] = W[2]
    c[t] = W[3]
    

    l[t] = (1/(ϕ + α)) * a[t] + (α/(ϕ + α)) * k[t] - (θ/(ϕ + α)) * c[t]
    

    y[t] = a[t] + α * k[t] + (1 - α) * l[t]
end

periods_to_report = [1, 4, 8, 16]
println("Output (y_t) responses:")
for p in periods_to_report
    println("t = $p: ", round(y[p], digits=6))
end

plot(1:T, y, title="Output (y_t) over Time", xlabel="Time (t)", ylabel="Output (y_t)", linewidth=2, color=:blue)