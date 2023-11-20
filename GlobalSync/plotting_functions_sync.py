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

    # Matrici JF e JH
    JF = np.array([[-2 * sigma_re, 0], [-2 * beta_im * sigma_re / beta_re, 0]])
    JH = np.array([[z0 ** (m - 1) * m * mu_re, -z0 ** (m - 1) * mu_im],
                   [z0 ** (m - 1) * m * mu_im, z0 ** (m - 1) * mu_re]])

    # Calcolo della legge di dispersione
    nx = 500
    xx = np.linspace(np.min(np.real(Lambdak)), 0, nx)
    cont_reldisp = np.zeros(nx)

    for ii in range(nx):
        xp = xx[ii]
        J = JF + xp * JH
        lambda_vals = np.linalg.eigvals(J)
        cont_reldisp[ii] = np.max(np.real(lambda_vals))

    # Calcolo per diversi Lambdak
    def calcola_dispersione(Lambda_values):
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


def plot_results(sol, k, config):
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

    tmin = 0
    tmax = ttime[-1]
    ixt = np.where((ttime > tmin) & (ttime < tmax))[0]

    # Prepare figure with multiple subplots
    fig, axs = plt.subplots(3, 1, dpi=200, figsize=(6, 8))
    # Set overall title
    fig.suptitle(f'Trajectories, k = {k}', fontfamily='serif', fontsize=16)
    # Plot real and imaginary parts of wt
    axs[0].plot(ttime[ixt], np.real(wt[:, ixt]).T, 'r')
    axs[0].plot(ttime[ixt], np.imag(wt[:, ixt]).T, 'b')
    axs[0].set_xlabel('time')
    axs[0].set_ylabel(r'$\mathcal{Re}, \mathcal{Im}$ of $w(t)$')

    line_re = mlines.Line2D([], [], color='red', label=r'$\mathcal{Re}$')
    line_im = mlines.Line2D([], [], color='blue', label=r'$\mathcal{Im}$')
    axs[0].legend(handles=[line_re, line_im], loc='upper right')


    # Heatmap of real part of wt
    im1 = axs[1].imshow(np.real(wt[:, ixt]), aspect='auto', extent=[ttime[ixt].min(), ttime[ixt].max(), 1, wt.shape[0]])
    fig.colorbar(im1, ax=axs[1])
    axs[1].set_xlabel('time')
    axs[1].set_ylabel(r'Real Part of $w(t)$')

    # Heatmap of magnitude of wt
    im2 = axs[2].imshow(np.abs(wt[:, ixt]), aspect='auto', extent=[ttime[ixt].min(), ttime[ixt].max(), 1, wt.shape[0]])
    fig.colorbar(im2, ax=axs[2])
    axs[2].set_xlabel('time')
    axs[2].set_ylabel('$|w(t)|$')

    plt.tight_layout()
    plt.show()

    num_subplots = 1  # Base per R e Rw
    if k != 0 and Bk is not None:
        num_subplots += 1
    if k != 3 and Bkp1 is not None:
        num_subplots += 1

    fig2, axs = plt.subplots(num_subplots, 1, dpi=200, figsize=(5, num_subplots * 2.8))

    if num_subplots == 1:
        axs = [axs]  # Assicurati che axs sia una lista

    # Compute R and Rw
    R = np.zeros(nt)
    Rw = np.zeros(nt)
    for ii in range(nt):
        R[ii] = np.abs(np.sum(np.exp(1j * np.angle(wt[:, ii]))) / wt.shape[0])
        Rw[ii] = np.abs(np.sum(wt[:, ii]) / wt.shape[0])

    # Plot R and Rw
    axs[0].plot(ttime, R, label='R(t)', linewidth=2)
    axs[0].plot(ttime, Rw, label='Rw(t)', linewidth=2)
    axs[0].set_xlabel('time')
    axs[0].set_ylabel('$R(t)$ and $R_w(t)$')
    axs[0].legend()
    axs[0].set_title(f'Order parameters, k = {k}', fontfamily='serif', fontsize=13)
    # Index per tracciare i subplot aggiuntivi
    subplot_index = 1

    # Additional plots for Bkwt and Bkp1wt if k is specified
    if k != 0 and Bk is not None:
        Bkwt = Bk @ wt
        axs[subplot_index].plot(ttime[ixt], np.real(Bkwt[:, ixt]).T, 'r')
        axs[subplot_index].plot(ttime[ixt], np.imag(Bkwt[:, ixt]).T, 'b')
        axs[subplot_index].set_xlabel('time')
        axs[subplot_index].set_ylabel(f'$B_{{{k}}}$' + '$w(t)$')
        line_re = mlines.Line2D([], [], color='red', label=r'$\mathcal{Re}$')
        line_im = mlines.Line2D([], [], color='blue', label=r'$\mathcal{Im}$')
        axs[subplot_index].legend(handles=[line_re, line_im], loc='upper right')
        axs[subplot_index].set_title(f'Projected down signal', fontfamily='serif', fontsize=13)
        subplot_index += 1

    if k != 3 and Bkp1 is not None:
        Bkp1wt = Bkp1.T @ wt
        axs[subplot_index].plot(ttime[ixt], np.real(Bkp1wt[:, ixt]).T, 'r')
        axs[subplot_index].plot(ttime[ixt], np.imag(Bkp1wt[:, ixt]).T, 'b')
        axs[subplot_index].set_xlabel('time')
        axs[subplot_index].set_ylabel(f'$B_{{{k+1}}}^T$' + '$w(t)$')
        line_re = mlines.Line2D([], [], color='red', label=r'$\mathcal{Re}$')
        line_im = mlines.Line2D([], [], color='blue', label=r'$\mathcal{Im}$')
        axs[subplot_index].legend(handles=[line_re, line_im], loc='upper right')
        axs[subplot_index].set_title(f'Projected up signal', fontfamily='serif', fontsize=13)

    plt.tight_layout()
    plt.show()
