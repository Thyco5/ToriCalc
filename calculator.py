# -*- coding: utf-8 -*-
import math

def calculate_chapman_enskog_DAB(T, P, M_A, M_B, sigma_A, sigma_B, epsilon_k_A, epsilon_k_B):
    """
    Calculate the binary gas diffusion coefficient (D_AB) at low pressures using the Chapman-Enskog theory (Eq. 11-3.2).

    Parameters:
        T (float): Temperature [K]
        P (float): Pressure [bar]
        M_A (float): Molar mass of component A [g/mol]
        M_B (float): Molar mass of component B [g/mol]
        sigma_A (float): Lennard-Jones size parameter for A [Å]
        sigma_B (float): Lennard-Jones size parameter for B [Å]
        epsilon_k_A (float): Lennard-Jones energy parameter for A (epsilon/k) [K]
        epsilon_k_B (float): Lennard-Jones energy parameter for B (epsilon/k) [K]

    Returns:
        float: D_AB [cm^2/s]
    """
    # Error handling
    if T <= 0:
        raise ValueError("Temperature must be positive.")
    if P <= 0:
        raise ValueError("Pressure must be positive.")
    if sigma_A <= 0 or sigma_B <= 0:
        raise ValueError("LJ sigma values must be positive.")
    if epsilon_k_A <= 0 or epsilon_k_B <= 0:
        raise ValueError("LJ epsilon/k values must be positive.")
    if M_A <= 0 or M_B <= 0:
        raise ValueError("Molar masses must be positive.")

    # Step 1: Calculate M_AB
    M_AB = 2 / (1/M_A + 1/M_B)

    # Step 2: Calculate sigma_AB
    sigma_AB = (sigma_A + sigma_B) / 2

    # Step 3: Calculate epsilon_AB/k
    epsilon_AB_k = math.sqrt(epsilon_k_A * epsilon_k_B)

    # Step 4: Calculate T_star
    T_star = T / epsilon_AB_k

    # Step 5: Calculate Omega_D (Neufeld et al. relation)
    A = 1.06036
    B = 0.15610
    C = 0.19300
    D = 0.47635
    E = 1.03587
    F = 1.52996
    G = 1.76474
    H = 3.89411
    Omega_D = (A / (T_star ** B)) \
        + (C / math.exp(D * T_star)) \
        + (E / math.exp(F * T_star)) \
        + (G / math.exp(H * T_star))

    # Step 6: Calculate D_AB (Eq. 11-3.2)
    D_AB = (0.00266 * T ** 1.5) / (P * (M_AB ** 0.5) * (sigma_AB ** 2) * Omega_D)

    return D_AB

# Example test (N2-CO2 at 590 K, 1 bar)
if __name__ == "__main__":
    T = 590  # K
    P = 1    # bar
    M_A = 28.0  # g/mol (N2)
    M_B = 44.0  # g/mol (CO2)
    sigma_A = 3.798  # Å (N2)
    sigma_B = 3.941  # Å (CO2)
    epsilon_k_A = 71.4   # K (N2)
    epsilon_k_B = 195.2  # K (CO2)
    D_AB = calculate_chapman_enskog_DAB(T, P, M_A, M_B, sigma_A, sigma_B, epsilon_k_A, epsilon_k_B)
    print("D_AB = {:.3f} cm^2/s (Expected approx. 0.52 cm^2/s)".format(D_AB))
