import numpy as np
from scipy.integrate import solve_ivp
from functions_turing import *
from plotting_functions_turing import plot_aggregated_results, calculate_and_plot_dispersion
"""
Configuration file for the Turing Pattern simulation. The one named 'Muolo_etal.yaml' is the one used
in the cross-diffusion case and is used in the paper `The three way Dirac operator and dynamical Turing and Dirac 
induced patterns on nodes and links` by Riccardo Muolo, Timoteo Carletti and Ginestra Bianconi. 
The one named 'Giambagli_etal_PRE.yaml' is the one used to reproduce the results of the paper `Diffusion-driven 
instability of topological signals coupled by the Dirac operator` by Lorenzo Giambagli, Lucille Calmon, Riccardo Muolo, 
Timoteo Carletti and Ginestra Bianconi.
"""

configuration_name = 'Muolo_etal.yaml'

# Load the configuration file
config = LoadConfig(configuration_name)
A, B = import_matrices('Matrices_benchmark.npz')['adjacency'], import_matrices('Matrices_benchmark.npz')['boundary']

# Define the hodge Laplacian and other initializations
nodes, links = B.shape

# Hodge Laplacian, nodes and links
deg = np.sum(A, axis=1)  # Sum of A along rows
L0 = np.diag(deg) - A    # L0 as the difference of diagonal degree matrix and adjacency matrix A
L1 = B.T @ B  # Multiplication of B transpose with B

# DISPERSION RELATION
calculate_and_plot_dispersion(L0, config)

# INTEGRATION
t_ini, t_final, y0 = initial_condition(nodes, links, config)

# Perform the integration using solve_ivp with the 'RK45' method, which is equivalent to RK4 for our purposes
sol = solve_ivp(system_to_integrate, [t_ini, t_final], y0, method='RK45', args=(config, B, L0, L1))

# PLOTTING
plot_aggregated_results(sol, B)


