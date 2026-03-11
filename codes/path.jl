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
    xlabels="Okres", 
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
savefig(p, "sciezka_ar1_tauchen.png")

println("Symulacja zakończona. Plik 'sciezka_ar1_tauchen.jpg' został zapisany.")
println("Wartości stanów: ", round.(mc.state_values, digits=3))

# 2. Dyskretyzacja metodą Rouwenhorsta
# Uwaga: rouwenhorst w QuantEcon.jl generuje proces o średniej 0
# Parametry: n, ρ, σ_ϵ
mc_rouwenhorst_raw = rouwenhorst(n_states, ρ, σ_ϵ)

# Przesuwamy stany o μ, aby proces startował i oscylował wokół 1.0
states_shifted = mc_rouwenhorst_raw.state_values .+ μ
mc_rouwenhorst = MarkovChain(mc_rouwenhorst_raw.p, states_shifted)

# 3. Symulacja ścieżki
start_index = findmin(abs.(mc_rouwenhorst.state_values .- 1.0))[2]

# Używamy tego samego seeda co wcześniej, by porównać metody na tym samym "szumie"
Random.seed!(42) 
path_indices = simulate_indices(mc_rouwenhorst, n_periods; init=start_index)
path_values = mc_rouwenhorst.state_values[path_indices]

# 4. Generowanie wykresu
p2 = plot(path_values, 
    title="Symulacja procesu AR(1) - Metoda Rouwenhorsta",
    xlabel="Okres", 
    ylabel="Wartość stanu",
    label="Ścieżka procesu (Rouwenhorst)",
    marker=:square,
    markersize=3,
    color=:red,
    grid=true,
    ylims=(μ - 4*σ_y, μ + 4*σ_y))

hline!(mc_rouwenhorst.state_values, color=:gray, linestyle=:dot, alpha=0.5, label="Stany")

# 5. Zapis do pliku
savefig(p2, "sciezka_ar1_rouwenhorst.png")