import matplotlib.pyplot as plt
import numpy as np


def plot_aggregated_results(sol, B):
    """
    This function plots the results of the integration and the projection in the image of B and B.T of them.
    :param sol: Result of the integration using the scipy.integrate.solve_ivp function
    :param B: The boundary operator of the simplicial complex
    :return: None
    """
    ttime = sol.t
    nodes, _ = B.shape
    ut = sol.y[:nodes, :]
    vt = sol.y[nodes:, :]

    # Configure plot to use LaTeX for rendering if desired
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times', size=24, style='italic')

    fig, axs = plt.subplots(3, 2, figsize=(16, 16))  # Adjust the figsize to fit your needs

    # Plot u(t) over time
    for i in range(ut.shape[0]):
        axs[0, 0].plot(ttime, ut[i,:], linewidth=2)
    axs[0, 0].set_xlabel('time')
    axs[0, 0].set_ylabel(r'$u_i(t)$')

    # Plot v(t) over time
    for i in range(vt.shape[0]):
        axs[0, 1].plot(ttime, vt[i,:], linewidth=2)
    axs[0, 1].set_xlabel('time')
    axs[0, 1].set_ylabel(r'$v_i(t)$')

    # # Plot B' * u(t) over time
    projected_ut = B.T @ ut
    for i in range(projected_ut.shape[0]):
        axs[1, 0].plot(ttime, projected_ut[i,:], linewidth=2)
    axs[1, 0].set_xlabel('time')
    axs[1, 0].set_ylabel(r'$(B^\top u)_i(t)$')

    # Plot B * v(t) over time
    projected_vt = B @ vt
    for i in range(projected_vt.shape[0]):
        axs[1, 1].plot(ttime, projected_vt[i,:], linewidth=2)
    axs[1, 1].set_xlabel('time')
    axs[1, 1].set_ylabel(r'$(Bv)_i(t)$')

    # Heatmap for u(t)
    im0 = axs[2, 0].imshow(ut, aspect='auto', origin='lower')
    fig.colorbar(im0, ax=axs[2, 0])
    axs[2, 0].set_xlabel('t')
    axs[2, 0].set_ylabel('nodes')
    axs[2, 0].set_title('u')

    # Heatmap for v(t)
    im1 = axs[2, 1].imshow(vt, aspect='auto', origin='lower')
    fig.colorbar(im1, ax=axs[2, 1])
    axs[2, 1].set_xlabel('t')
    axs[2, 1].set_ylabel('links')
    axs[2, 1].set_title('v')

    plt.tight_layout()
    plt.show()

# Example usage:
# plot_aggregated_results(sol, B)


def calculate_and_plot_dispersion(L0, config, step=0.001):
    """
    Calculate the dispersion relation for the dynamical system of function.system_of_equations.
    All the paraemters are extracted from the configuration file apart from L0 and are the dynamical parameters of the system.
    L0 is the Laplacian of the simplicial complex which is isospectral to L1, in this case.
    :param step: Minimal step for the calculation of hte continuous dispersion law. Default value is 0.001.
    Smaller is more precise. THIS FUNCTION IS MODEL DEPENDENT.
    """

    # Define and parse the model_parameters
    a = config['model_parameters']['a']
    alpha = config['model_parameters']['alpha']
    b = config['model_parameters']['b']         # Useless due to the fact that the fixed point is 0
    beta = config['model_parameters']['beta']   # Useless due to the fact that the fixed point is 0
    c = config['model_parameters']['c']
    gamma = config['model_parameters']['gamma']

    # Define and parse the diffusion_coefficients
    D0 = config['diffusion_coefficients'].get('D0', 0)
    D1 = config['diffusion_coefficients'].get('D1', 0)
    D01 = config['diffusion_coefficients'].get('D01', 0)
    D10 = config['diffusion_coefficients'].get('D10', 0)

    # Calculate the eigenvalues of the Laplacian L0
    L_eig = np.linalg.eigvals(L0)
    UP = np.max(np.abs(L_eig))

    # Initialize arrays to store the dispersion relation data
    K = []
    R = []
    J = np.array([[-a, 0],
                  [0, -alpha]])

    # First loop to populate K and R
    k_values = np.arange(0, np.sqrt(UP) + 1, step)
    for k in k_values:
        D_k = np.array([[-k**2 * D0, k * c - k**3 * D01],
                        [k * gamma - k**3 * D10, -k**2 * D1]])

        J_tilde = J + D_k
        eigenvalues = np.linalg.eigvals(J_tilde)
        K.append(k**2)
        R.append(np.max(np.real(eigenvalues)))

    # Convert to numpy arrays for indexing in plotting
    K = np.array(K)
    R = np.array(R)

    # Second loop to calculate the dispersion based on L_eig
    R_eig = []
    for j in range(len(L_eig)):
        if abs(L_eig[j]) < 1E-10:
            L_eig[j] = 0
        D_eig = np.array([[-L_eig[j] * D0, np.sqrt(L_eig[j]) * c - L_eig[j]**(3/2) * D01],
                          [np.sqrt(L_eig[j]) * gamma - L_eig[j]**(3/2) * D10, -L_eig[j] * D1]])

        J_tilde_eig = J + D_eig
        lam = np.linalg.eigvals(J_tilde_eig)
        R_eig.append(np.max(np.real(lam)))

    # Plot the results
    plt.figure()
    plt.plot(K, R, '-b', linewidth=2)
    plt.scatter(np.real(L_eig), R_eig, color='r', s=100, facecolor='c', edgecolor='r')
    plt.plot(np.arange(0, int(np.max(K)) + 1), np.zeros(int(np.max(K)) + 1), '-k')
    plt.xlabel('$\Lambda_0$', fontsize=18)
    plt.ylabel('$\Re\lambda$', fontsize=18)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('Dispersion Relation', fontsize=20)
    plt.legend(['R(k)', 'R_eig'])
    plt.grid(True)
    plt.tight_layout()
    plt.show()

