using LinearAlgebra
using Plots
using Interpolations

# 1. Definicja parametrów (zgodnie z poprzednim zadaniem)
P       = [0.5 0.5; 0.5 0.5]
Y       = [0.4, 1.0]
M       = length(Y)
phi     = 0.0
N       = 100
Amax    = 20.0
A       = collect(range(-phi, Amax, length=N)) # Egzogenny grid dla a'
sigma   = 1.0
beta    = 0.99
R       = 0.99 * (1/beta)

# Pochodna funkcji użyteczności u'(c) i jej odwrotność (u')^-1
u_prime(c) = c^(-sigma)
u_prime_inv(val) = val^(-1/sigma)

# 2. Inicjalizacja
# Zaczynamy od zgadnięcia konsumpcji c(a, y) = R*a + y [cite: 116]
c_old = zeros(N, M)
for j in 1:M, i in 1:N
    c_old[i, j] = R * A[i] + Y[j]
end

c_new = copy(c_old)
tol = 1e-8
dist = 1.0
iter = 0

# 3. Pętla EGM
while dist > tol
    global iter += 1
    
    # KROK 1: Obliczenie oczekiwanej krańcowej użyteczności [cite: 118]
    # RHS = beta * R * E[u'(c(a', y'))]
    rhs = beta * R * (u_prime.(c_old) * P')
    
    # KROK 2: Wyznaczenie "quasi-polityki" konsumpcji c_bar [cite: 118, 119]
    c_bar = u_prime_inv.(rhs)
    
    for j in 1:M
        # KROK 3: Wyznaczenie endogennego gridu aktywów a_bar 
        # a_bar = (c_bar + a' - y) / R
        a_bar = (c_bar[:, j] .+ A .- Y[j]) ./ R
        
        # KROK 4: Mapowanie na stały grid egzogenny A [cite: 124]
        # (a) Obsługa ograniczenia kredytowego (a < a_bar_1) [cite: 125, 126, 127]
        # Jeśli aktywa są mniejsze niż pierwszy punkt endogennego gridu, agent pożycza do limitu a' = -phi
        
        # (c) Interpolacja liniowa dla punktów wewnątrz gridu [cite: 130, 131, 132]
        interp_func = LinearInterpolation(a_bar, c_bar[:, j], extrapolation_bc=Line())
        
        for i in 1:N
            if A[i] <= a_bar[1]
                # Przypadek (a): Ograniczenie wiąże [cite: 126]
                c_new[i, j] = R * A[i] + Y[j] - A[1]
            else
                # Przypadek (c): Interpolacja [cite: 132]
                c_new[i, j] = interp_func(A[i])
            end
        end
    end
    
    global dist = maximum(abs.(c_new .- c_old))
    c_old .= c_new
end

println("EGM zbieżny po $iter iteracjach.")

# 4. Generowanie wykresów
p1 = plot(A, c_old[:, 1], label="Niska prod. (Y=0.4)", title="Funkcja Konsumpcji (EGM)",
          xlabel="Aktywa (a)", ylabel="c(a, y)", lw=2)
plot!(A, c_old[:, 2], label="Wysoka prod. (Y=1.0)", lw=2)

# Obliczenie funkcji polityki aktywów a'
a_prime = [R * A[i] + Y[j] - c_old[i, j] for i in 1:N, j in 1:M]

p2 = plot(A, a_prime[:, 1], label="Niska prod. (Y=0.4)", title="Funkcja Polityki Aktywów a'(a, y)",
          xlabel="Aktywa (a)", ylabel="a'", lw=2)
plot!(A, a_prime[:, 2], label="Wysoka prod. (Y=1.0)", lw=2)
plot!(A, A, label="45 stopni", linestyle=:dash, color=:black)

savefig(p1, "egm_konsumpcja.png")
savefig(p2, "egm_aktywa.png")