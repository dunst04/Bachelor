using LinearAlgebra

# 1. Definicja parametrów
P       = [0.5 0.5; 0.5 0.5]      # Macierz przejścia dochodu [cite: 8, 87]
Y       = [0.4, 1.0]               # Grid dochodu 
M       = length(Y)
phi     = 0.0                      # Limit zadłużenia [cite: 85, 93]
N       = 100                      # Liczba punktów gridu aktywów
Amax    = 20.0
A       = collect(range(-phi, Amax, length=N)) # Grid aktywów 
sigma   = 1.0                      # Parametr użyteczności (CRRA: u(c) = c^(1-σ)/(1-σ))
beta    = 0.99                     # Czynnik dyskontujący [cite: 86]
R       = 0.99 * (1/beta)          # Stopa zwrotu [cite: 86]

# Funkcja użyteczności (CRRA)
function utility(c, σ)
    if c <= 0
        return -1e10 # Kara za ujemną konsumpcję
    end
    return σ == 1 ? log(c) : (c^(1 - σ)) / (1 - σ)
end

# 2. Inicjalizacja
V_old = zeros(N, M) # Początkowe zgadnięcie funkcji wartości [cite: 94]
V_new = zeros(N, M)
policy_a = zeros(Int, N, M) # Indeksy optymalnych aktywów a' [cite: 88]

tol = 1e-6       # Kryterium zbieżności [cite: 98]
dist = 1.0
iter = 0

# 3. Pętla VFI [cite: 95, 103]
while dist > tol
    global iter += 1
    
    # Obliczamy oczekiwaną wartość przyszłą dla każdego możliwego a' [cite: 84]
    # EV(a', y) = Σ P(y'|y) * V(a', y')
    EV = V_old * P' 

    for j in 1:M      # Pętla po obecnym dochodzie y
        for i in 1:N  # Pętla po obecnych aktywach a
            max_val = -1e12
            best_a_idx = 1
            
            for k in 1:N # Pętla po wyborze aktywów na przyszły okres a'
                # Ograniczenie budżetowe: c + a' = R*a + y [cite: 85, 96]
                cons = R * A[i] + Y[j] - A[k]
                
                if cons > 0
                    # Bellman: u(c) + β * E[V(a', y')] [cite: 84, 95]
                    val = utility(cons, sigma) + beta * EV[k, j]
                    
                    if val > max_val
                        max_val = val
                        best_a_idx = k
                    end
                end
            end
            V_new[i, j] = max_val
            policy_a[i, j] = best_a_idx
        end
    end

    # Sprawdzenie zbieżności (norma max) [cite: 98, 100]
    global dist = maximum(abs.(V_new .- V_old))
    V_old .= V_new
    
    if iter % 50 == 0
        println("Iteracja: $iter, Dystans: $dist")
    end
end

println("Zbieżność osiągnięta po $iter iteracjach.")

using Plots

# --- [Uruchomienie obliczeń VFI z poprzedniego kroku] ---
# Zakładamy, że zmienne V_old, policy_a, A, Y, R i sigma są już w pamięci

# 4. Obliczenie funkcji konsumpcji
# c(a, y) = R*a + y - a'(a, y)
consumption = zeros(N, M)
for j in 1:M
    for i in 1:N
        a_prime = A[policy_a[i, j]]
        consumption[i, j] = R * A[i] + Y[j] - a_prime
    end
end

# 5. Generowanie wykresu Funkcji Wartości
p1 = plot(A, V_old[:, 1], 
    label="Niska prod. (Y=0.4)", 
    title="Funkcja Wartości V(a, y)",
    xlabel="Aktywa (a)", 
    ylabel="Wartość V",
    linewidth=2, color=:blue)
plot!(A, V_old[:, 2], 
    label="Wysoka prod. (Y=1.0)", 
    linewidth=2, color=:red)

# 6. Generowanie wykresu Funkcji Konsumpcji
p2 = plot(A, consumption[:, 1], 
    label="Niska prod. (Y=0.4)", 
    title="Funkcja Konsumpcji c(a, y)",
    xlabel="Aktywa (a)", 
    ylabel="Konsumpcja c",
    linewidth=2, color=:blue,
    legend=:topleft)
plot!(A, consumption[:, 2], 
    label="Wysoka prod. (Y=1.0)", 
    linewidth=2, color=:red)

# 7. Zapis do plików PNG
savefig(p1, "funkcja_wartosci.png")
savefig(p2, "funkcja_konsumpcji.png")

println("Wykresy zostały zapisane jako funkcja_wartosci.png i funkcja_konsumpcji.png")