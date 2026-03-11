using QuantEcon
using Plots
using Random

# 1. Definicja parametrów
ρ = 0.9                # Autokorelacja
σ_y = 0.1              # Bezwarunkowe odchylenie standardowe procesu
μ = 1.0                # Średnia bezwarunkowa (i punkt startowy)
n_states = 7           # Liczba stanów
n_periods = 100        # Czas trwania symulacji

# Obliczamy odchylenie standardowe szoku (white noise)
# Ponieważ Var(y) = Var(ϵ) / (1 - ρ^2)
σ_ϵ = σ_y * sqrt(1 - ρ^2)

# 2. Dyskretyzacja Tauchena
# Parametry: (liczba stanów, ρ, σ_ϵ, μ, n_std)
# n_std to liczba odchyleń standardowych określająca zakres siatki (standardowo 3)
mc = tauchen(n_states, ρ, σ_ϵ, (1 - ρ) * μ, 3)

# 3. Symulacja ścieżki
# Znajdujemy indeks stanu najbliższego wartości startowej 1.0
start_index = findmin(abs.(mc.state_values .- 1.0))[2]

Random.seed!(42) # Dla powtarzalności wyników
path_indices = simulate_indices(mc, n_periods; init=start_index)
path_values = mc.state_values[path_indices]

# 4. Generowanie wykresu
p = plot(path_values, 
    title="Symulacja procesu AR(1) - Algorytm Tauchena",
    xlabel="Okres", 
    ylabel="Wartość stanu",
    label="Ścieżka procesu",
    marker=:circle,
    markersize=3,
    color=:blue,
    grid=true,
    ylims=(μ - 3*σ_y, μ + 3*σ_y))

# Dodanie linii siatki dla dyskretnych stanów
hline!(mc.state_values, color=:gray, linestyle=:dash, alpha=0.5, label="Stany")

# 5. Zapis do pliku JPG
savefig(p, "sciezka_ar1_tauchen.jpg")

println("Symulacja zakończona. Plik 'sciezka_ar1_tauchen.jpg' został zapisany.")
println("Wartości stanów: ", round.(mc.state_values, digits=3))