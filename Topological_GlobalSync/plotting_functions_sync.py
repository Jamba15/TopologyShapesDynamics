import matplotlib.pyplot as plt
import numpy as np

from functions_sync import *
import matplotlib.lines as mlines


def calculate_and_plot_dispersion(config, k):

    Lk, Lkup, Lkdwn = get_Laplacian(k, config['simplex_boundary_filename'])

    # Eigenvalues and Eigenvectors
    Lambdak, VVk = np.linalg.eig(Lk)

    Lambdakup, VVkup = np.linalg.eig(Lkup)

    Lambdakdwn, VVkdwn = np.linalg.eig(Lkdwn)

    sigma = config['model_parameters']['sigma']
    beta = config['model_parameters']['beta']
    mu = config['model_parameters']['mu']
    m = config['model_parameters']['m']

    sigma_re = np.real(sigma)
    sigma_im = np.imag(sigma)
    beta_re = np.real(beta)
    beta_im = np.imag(beta)
    mu_re = np.real(mu)
    mu_im = np.imag(mu)

    z0 = np.sqrt(sigma_re / beta_re)

    # Calculate the Jacobian
    JF = np.array([[-2 * sigma_re, 0], [-2 * beta_im * sigma_re / beta_re, 0]])
    JH = np.array([[z0 ** (m - 1) * m * mu_re, -z0 ** (m - 1) * mu_im],
                   [z0 ** (m - 1) * m * mu_im, z0 ** (m - 1) * mu_re]])

    # Evaluate the dispersion relation
    nx = 500
    xx = np.linspace(np.min(np.real(Lambdak)), 0, nx)
    cont_reldisp = np.zeros(nx)

    for ii in range(nx):
        xp = xx[ii]
        J = JF + xp * JH
        lambda_vals = np.linalg.eigvals(J)
        cont_reldisp[ii] = np.max(np.real(lambda_vals))


    def calcola_dispersione(Lambda_values):
        """
        Evaluate the dispersion relation for a given set of eigenvalues
        :param Lambda_values: The eigenvalues to evaluate
        :return: The dispersion relation evaluated at the given eigenvalues
        """
        nk = Lambda_values.shape[0]
        dispersion_values = np.zeros(nk)
        for ii in range(nk):
            J = JF + Lambda_values[ii] * JH
            lambda_vals = np.linalg.eigvals(J)
            dispersion_values[ii] = np.max(np.real(lambda_vals))
        return dispersion_values

    simplex_reldisp_k = calcola_dispersione(Lambdak)
    simplex_reldisp_kup = calcola_dispersione(Lambdakup)
    simplex_reldisp_kdwn = calcola_dispersione(Lambdakdwn)

    # Creazione del grafico
    plt.plot(-xx, cont_reldisp, 'b-', linewidth=2)
    plt.plot(-xx, np.zeros_like(xx), 'k-', linewidth=2)

    plt.plot(-Lambdak, simplex_reldisp_k, 'ro', markerfacecolor='r',  alpha=0.5, markersize=8, label='$\Lambda_k$')
    if k != 0:
        plt.plot(-Lambdakup, simplex_reldisp_kup, 'gs', markerfacecolor='g',  alpha=0.5,  markersize=8, label='$\Lambda^{up}_{k}$')
    if k != 3:
        plt.plot(-Lambdakdwn, simplex_reldisp_kdwn, 'cd', markerfacecolor='c',  alpha=0.5,  markersize=8, label='$\Lambda^{down}_{k}$')

    plt.title('Dispersion relation', fontsize=24)
    plt.xlabel('r$-\Lambda^{(\\alpha)}$', fontsize=24)
    plt.ylabel('$MSF$', fontsize=24)
    plt.legend(fontsize=24)
    plt.tight_layout()
    plt.show()


def plot_trajectories(sol, k, config, **kwargs):
    # Extract B matrices according to k
    Bk = None
    Bkp1 = None
    if k != 0:
        Bk = import_matrices(config['simplex_boundary_filename'])[f'B{k}']
    if k != 3:
        Bkp1 = import_matrices(config['simplex_boundary_filename'])[f'B{k + 1}']

    ttime = sol.t
    wt = sol.y
    nt = len(ttime)

    tmin = kwargs.get('tmin', 50)
    tmax = ttime[-1]
    ixt = np.where((ttime > tmin) & (ttime < tmax))[0]

    # Prepare figure with multiple subplots
    fig, axs = plt.subplots(2, 1, dpi=200, figsize=(6, 8))
    # Set overall title
    fig.suptitle(f'Trajectories', fontfamily='serif', fontsize=16)

    # Plot real and imaginary parts of wt as two different plots
    axs[0].plot(ttime[ixt], np.real(wt[:, ixt]).T, 'r')
    axs[0].set_xlabel('time')
    axs[0].set_ylabel(r'$\mathcal{Re}(w(t))$')

    axs[1].plot(ttime[ixt], np.imag(wt[:, ixt]).T, 'b')
    axs[1].set_xlabel('time')
    axs[1].set_ylabel(r'$\mathcal{Im}(w(t))$')

    plt.tight_layout()
    plt.show()


def plot_order_parameter(sol):
    ttime = sol.t
    wt = sol.y
    nt = len(ttime)

    # Compute R and Rw
    R = np.zeros(nt, dtype=complex)
    Rw = np.zeros(nt)

    for ii in range(nt):
        R[ii] = (np.sum(np.abs(wt[:, ii]) * np.exp(1j * np.angle(wt[:, ii]))) / wt.shape[0])
        Rw[ii] = np.abs(np.sum(wt[:, ii]) / wt.shape[0])

    tsave = int(np.floor(nt * 4 / 5))
    # print(tsave, nt)

    fig, ax = plt.subplots(1, 2, dpi=300, figsize=(6, 3))

    # Plot X1

    ax[0].plot(np.abs(R), linewidth=2)
    ax[0].set_xlabel('$t$')
    ax[0].set_ylabel('$R=|X|$')
    ax[0].set_title(f'Real order parameter', fontfamily='serif', fontsize=13)

    ax[1].plot(np.real(R[-tsave:-1]), np.imag(R[-tsave:-1]), linewidth=2)

    ax[1].set_xlabel('$Re(X)$')
    ax[1].set_ylabel(' $Im(X)$')
    ax[1].set_title(f'Complex order parameter', fontfamily='serif', fontsize=13)

    plt.tight_layout()
    plt.show()
    
    
    