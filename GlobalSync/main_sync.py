from functions_sync import *
from plotting_functions_sync import calculate_and_plot_dispersion, plot_results
from scipy.integrate import solve_ivp

config = LoadConfig('PRL.yml')

# Define the level of the simplicial complex. k = 0 is the vertex level, k = 1 is the edge level, etc.
# In the present version of the code k < 4
k = 3

# DISPERSION REALATION
calculate_and_plot_dispersion(config, k)

# INTEGRATION
t_ini, t_final, w_init, Lk = initial_condition(k, config)

# Integrate the system of differential equations with the Runge-Kutta 4 method
sol = solve_ivp(system_to_integrate, [t_ini, t_final], w_init, method='RK45', args=(config, Lk))

# PLOTTING
plot_results(sol, k, config)
