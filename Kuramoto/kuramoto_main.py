from functions_kuramoto import *
Number_of_filled_cells = 1
s_c = 1 # concordant

# Create octogons
x,y = octogon_positions()
a, OmegaE = octogon_adjacency(s_c)
I, J = np.where(np.triu(a))
I2, J2 = np.where(np.triu(OmegaE))
L = len(I) #number of edges
# Define internal frequencies:
Omega = (5*OmegaE[I2, J2] + 3*np.random.randn(L)).T

# Boundary matrices:
B1 = create_B1(L,I,J)
B2, B2b = create_B2(L,I,J, Number_of_filled_cells)

plot_octogons(x, y, I, J)
plt.show()

#%%
config = LoadConfig("kuramoto_config.yml")
config['Number_of_filled_cells'] = Number_of_filled_cells
config['s_c'] = s_c

def kuramoto_integrate_and_plot(Psi0, config, B1, B2, B2b, Omega, create_animation=True):
    Tmax = config['Tmax']  # max time to integrate
    t_save = config['t_save']  # transient time
    dt = config['dt']  # time step
    sigma = config['sigma']  #coupling strength

    x, y = octogon_positions()
    a, OmegaE = octogon_adjacency(s_c)
    I, J = np.where(np.triu(a))

    # Auxiliary variables for plots
    it_save = int(t_save / dt)
    tt = np.arange(0, Tmax, dt)  # array of times
    nts = len(tt)  # number of time-steps
    nts_save = nts - it_save  # number of time-steps to save

    # Order parameters
    R1a = np.zeros(nts_save, dtype=complex)  # synchronization of first octogon
    R2a = np.zeros(nts_save, dtype=complex)  # synchronization of second octogon
    idx_1 = np.asarray([0, 1, 3, 5, 6, 7, 8, 9])  # indices of edges in first octogon
    idx_2 = np.asarray([0, 2, 4, 10, 11, 12, 13, 14])  # indices of edges in second octogon

    # ---  Simulation ---
    Psi = Psi0.copy()  # Create a copy to modify

    for it in tqdm(range(nts), disable=not create_animation):  # Optional progress bar
        Psi = kuramoto_update(Psi, dt, sigma, B1, B2, Omega)
        Psi = np.mod(Psi, 2 * np.pi)

        if it >= it_save:
            x1 = B2b[idx_1, 0] * Psi[idx_1]
            R1a[it - it_save] = np.sum(np.exp(1j * x1)) / 8
            x2 = B2b[idx_2, 1] * Psi[idx_2]
            R2a[it - it_save] = np.sum(np.exp(-1j * x2)) / 8

        if create_animation and it % 20 == 0:
            fig = plt.figure(figsize=(6, 3), dpi=100)
            plot_octogons(x, y, I, J, Psi, save=True, filename=f'octogons_{it}')
            plt.close(fig)

    if create_animation:
        this_dir = os.path.dirname(os.path.abspath(__file__))
        input_directory = join(this_dir, 'animations')
        output_gif_name = "kuramoto_octagons.gif"
        create_gif(input_directory, output_gif_name)
    return R1a, R2a

initial_conditions = np.random.rand(B1.shape[1]) * 2 * np.pi
kuramoto_integrate_and_plot(initial_conditions, config, B1, B2, Omega, create_animation=True)
