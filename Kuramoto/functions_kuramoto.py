##%% Import modules and define functions
import os

import imageio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.colors as colors
import fnmatch
import os
from os.path import join, dirname, abspath
import yaml
import numpy as np
import re

from tqdm import tqdm


def Kuramoto_edges(Phi, sigma, B1, B2, Omega):
    """
    This function calculates the time derivative of the edge phases for the simplicial Kuramoto model  defined on the edges
    of a cell complex. As inputs it requires:
    Phi: previous state of the system, numpy array of dimension E = number of edges
    sigma: coupling strength (escalar)
    B1: boundary matrix of the cell complex at dimension 1, numpy array of dimension N x E, where N is the number of nodes.
    B2: boundary matrix of the cell complex at dimension 2, numpy array of dimension E x Nt, where Nt is the number of 2D cells
        or polygons.
    Omega: edge intrinsic frequency, numpy array of dimension E.
    """
    if len(B2.shape) > 1:
        Phi_dot = Omega - sigma * (
                np.matmul(B1.T, np.sin(np.matmul(B1, Phi)))
                + np.matmul(B2, np.sin(np.matmul(B2.T, Phi))))
    else:
        Phi_dot = Omega - sigma * (np.matmul(B1.T, np.sin(np.matmul(B1, Phi))) +
                                   B2 * np.sin(B2 * Phi))
    return Phi_dot


def octogon_positions():
    """
    Auxiliary function to define the location of the nodes (x,y) of the cell complex made by two adjacent octogons.
    For illustration purposes only.
    """

    N = 14
    x = np.zeros(N)
    y = np.zeros(N)
    x[:8] = np.cos(2 * np.pi * np.arange(1, 9) / 8)
    y[:8] = np.sin(2 * np.pi * np.arange(1, 9) / 8)
    x0 = x[0] + x[1]
    y0 = y[0] + y[1]
    theta1 = np.arctan((y[0] - y0) / (x[0] - x0))
    x[8:] = x0 + np.cos(-theta1 + 2 * np.pi * np.arange(1, 7) / 8)
    y[8:] = y0 + np.sin(-theta1 + 2 * np.pi * np.arange(1, 7) / 8)

    return x, y


def octogon_adjacency(s_c):
    """
    Function to define the adjacency matrix A (dimension: NxN) of the double-octogon cell complex.
    The auxiliary matrix OmegaE (dimension: NxN) which encodes the relative sign of the intrinsic frequencies of the edges
    is also defined. This matrix will be useful to define the vector of intrinsic frequencies.
    """
    N = 14
    A = np.zeros((N, N))
    OmegaE = np.zeros((N, N))

    for i in range(7):
        A[i, i + 1] = 1
        A[i + 1, i] = 1

        OmegaE[i, i + 1] = -s_c
        OmegaE[i + 1, i] = -s_c

    A[7, 0] = 1
    A[0, 7] = 1
    OmegaE[7, 0] = s_c;
    OmegaE[0, 7] = s_c;

    for i in range(8, 13):
        A[i, i + 1] = 1
        A[i + 1, i] = 1
        OmegaE[i, i + 1] = 1;
        OmegaE[i + 1, i] = 1;

    A[0, 8] = 1
    A[8, 0] = 1
    OmegaE[0, 8] = 1
    OmegaE[8, 0] = 1

    A[13, 1] = 1
    A[1, 13] = 1
    OmegaE[13, 1] = -1
    OmegaE[1, 13] = -1

    return A, OmegaE


def create_B1(L, I, J):
    """
    Funcion to define the B1 boundary matrix of the double-octogon cell complex. The inputs are:
    L: number of edges ( = E)
    I, J: indices of existing edges in the adjacency matrix.
    """
    N = 14
    B1 = np.zeros((N, L))

    for n in range(L):
        B1[I[n], n] = -1
        B1[J[n], n] = 1

    return B1


def create_B2(L, I, J, Number_of_filled_cells):
    """
    Funcion to define the B2 boundary matrix of the double-octogon cell complex. The inputs are:
    L: number of edges ( = E)
    I, J: indices of existing edges in the adjacency matrix
    Number_of_filled_cells: as the name indicates, number of octogons that are filled. Options are 0, 1, 2

    The outputs are B2 and the auxiliary matrix B2b, which is the boundary of the second octogon,
    and is used in the calculation of the order parameters.
    """

    B2a = np.zeros((L, 2))  # Boundary of first octogon
    B2b = np.zeros((L, 2))  # Boundary of second octogon

    for inn in range(7):
        n1 = np.where((I == inn) * (J == (inn + 1)))[0][0]
        B2a[n1, 0] = 1
    n1 = np.where((I == 0) * (J == 7))[0][0]
    B2a[n1, 0] = -1

    B2b[:, 0] = B2a[:, 0]
    for inn in range(8, 13):
        n1 = np.where((I == inn) * (J == (inn + 1)))[0][0]
        B2b[n1, 1] = 1
    n1 = np.where((I == 0) * (J == 8))[0][0]
    B2b[n1, 1] = 1
    n1 = np.where((I == 1) * (J == 13))[0][0]
    B2b[n1, 1] = -1
    n1 = np.where((I == 0) * (J == 1))[0][0]
    B2b[n1, 1] = -1

    if Number_of_filled_cells == 1:
        B2 = B2a
    elif Number_of_filled_cells == 2:
        B2 = B2b
    else:
        B2 = np.asarray([0])
    return B2, B2b

def LoadConfig(config_name, hyper=None):
    """
    This function loads the configuration file from the configuration folder. If hyper is not None, it will replace the
    hyperparameters in the configuration file with the ones in hyper
    :param config_name: The name of the configuration file
    :param hyper: The dictionary with the hyperparameters to replace in the configuration file
    :return: The configuration dictionary
    """
    # Check if hyper is a dictionary of length 1, otherwise raise an error
    if hyper is not None:
        if not isinstance(hyper, dict):
            raise TypeError("Hyperparameters must be a dictionary")
        elif len(hyper) != 1:
            raise ValueError("Hyperparameters must be a dictionary of length 1")
    try:
        configuration_file = find(config_name, join(dirname(abspath(__file__)), 'configurations'))[0]
    except IndexError:
        raise FileNotFoundError(
            "Configuration file not found, check the name. Files in the configurations folder are:\n"
            " {}".format(os.listdir(join(dirname(abspath(__file__)), 'configurations'))))

    with open(configuration_file, 'r') as c:
        configuration = yaml.load(c, Loader=PrettySafeLoader)

    if hyper is not None:
        if type(hyper) is list:
            # if hyper value is of length larger than 1 pick the first one
            key, values = list(hyper.items())[0]
            if len(values) > 1:
                new_hyper = {key: values[0]}
            else:
                new_hyper = hyper
        else:
            new_hyper = hyper
        return DictionaryReplace(configuration, new_hyper)
    else:
        return configuration

class PrettySafeLoader(yaml.SafeLoader):
    def construct_python_tuple(self, node):
        return tuple(self.construct_sequence(node))

    # create constructor for !path tag
    def construct_path(self, node):
        return os.path.normpath(self.construct_scalar(node))

    def construct_complex(self, node):
        value = self.construct_scalar(node)
        # Utilizzare espressioni regolari per estrarre le parti reale e immaginaria
        match = re.match(r"([+-]?[0-9.]+)\s*([+-])\s*i([0-9.]+)", value)
        if match:
            real, sign, imag = match.groups()
            return complex(float(real), float(imag) if sign == '+' else -float(imag))
        else:
            raise ValueError(f"Non è possibile convertire '{value}' in un numero complesso")


def DictionaryReplace(dictionary: dict, sobstitute: dict) -> dict:
    """
    Sostituisce le chiavi di un dizionario con i valori di un altro dizionario. La sostituizione avviene in profondità
    e solo se la chiave è uguale
    :param dictionary: Il primo dizionario in cui la sostituzione deve avvenire in accordo alle chiavi di sobstitute
    :param sobstitute: Dizionario con le chiavi da sostituire
    :return: Il dizionario con i valori sostituiti in accordo alle chiavi di sobstitute
    """
    for key, value in sobstitute.items():
        for original_key, original_value in dictionary.items():
            if isinstance(original_value, dict):
                dictionary[original_key] = DictionaryReplace(original_value, sobstitute)
            elif original_key == key:
                dictionary[original_key] = value
    return dictionary

def kuramoto_update(Psi, dt, sigma, B1, B2, Omega):
    """Performs one RK4 update step for the Kuramoto model."""
    k1 = Kuramoto_edges(Psi, sigma, B1, B2, Omega)
    k2 = Kuramoto_edges(Psi + 0.5 * dt * k1, sigma, B1, B2, Omega)
    k3 = Kuramoto_edges(Psi + 0.5 * dt * k2, sigma, B1, B2, Omega)
    k4 = Kuramoto_edges(Psi + dt * k3, sigma, B1, B2, Omega)
    return Psi + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

PrettySafeLoader.add_constructor(
    u'tag:yaml.org,2002:python/tuple',
    PrettySafeLoader.construct_python_tuple)

PrettySafeLoader.add_constructor(
    '!path',
    PrettySafeLoader.construct_path)

PrettySafeLoader.add_constructor(
    '!complex',
    PrettySafeLoader.construct_complex)


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):

        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def plot_octogons(x, y, I, J, Psi=None, colormap="rainbow", save=False, filename="octogons.png"):
    """Plots the octogon structure with edges color-coded according to their phase.

    Args:
        x, y: tuple of numpy.ndarray (Octogon coordinates)
        I, J: numpy.ndarray (Edge indices)
        Psi:  numpy.ndarray (Phase values of the edges)
        cc:   str or matplotlib.colors.Colormap (Colormap for edges)
    """

    cc_off = (255. / 255., 171. / 255., 91. / 255., 1)  # Node color

    if Psi is not None:
        plt.title("Phase of the edges")
        # Normalize Psi for colormapping
        Psi_normalized = (Psi + np.pi) / (2 * np.pi)  # Map [-pi, pi] to [0, 1]
    else:
        # Set title to "Plotting the network skeleton of the cell complex" 
        plt.title("The network skeleton of the cell complex")
        Psi_normalized = None

    # Check if cc is a string or a Colormap object
    if isinstance(colormap, str):
        cc = cm.get_cmap(colormap)  # Get Colormap object from name
    else:
        raise ValueError("cc must be a string")

    if Psi is not None:
        # Plot the edges with colors from the colormap
        for n in range(len(I)):
            color = cc(Psi_normalized[n])
            plt.plot([x[I[n]], x[J[n]]], [y[I[n]], y[J[n]]], color=color, linewidth=4)
    else:
        for n in range(len(I)):
            plt.plot([x[I[n]], x[J[n]]], [y[I[n]], y[J[n]]], color='black', linewidth=4)

    # Plot nodes with fixed color
    plt.scatter(x, y, color=cc_off, s=100, label='Nodes')
    plt.axis('off') 

    if save:
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(join(this_dir, 'animations')):
            os.makedirs(join(this_dir, 'animations'))
        # save in that folder with an incremental filename
        plt.savefig(join(this_dir, 'animations', filename))
        plt.close()


def create_gif(input_dir, gif_name):
    """Creates a GIF from PNG files in a directory.

    Args:
        input_dir (str): Path to the directory containing PNG images.
        gif_name (str): Name of the output GIF file.
    """

    images = []
    filenames = sorted([f for f in os.listdir(input_dir) if f.endswith('.png')])
    for filename in tqdm(filenames):
        images.append(imageio.v3.imread(os.path.join(input_dir, filename)))  # Using imread_v3
    # save in one folder above input_dir
    imageio.mimsave(os.path.join(input_dir, '..', gif_name), images)
    print('GIF saved as', os.path.join(input_dir, '..', gif_name))
    # remove PNG files
    for filename in filenames:
        os.remove(os.path.join(input_dir, filename))
    print('PNG files removed')


def kuramoto_integrate_and_plot(Psi0, config, B1, B2, B2b, Omega, create_animation=True):
    Tmax = config['Tmax']  # max time to integrate
    t_save = config['t_save']  # transient time
    dt = config['dt']  # time step
    sigma = config['sigma']  #coupling strength

    # Auxiliary variables for plots
    it_save = int(t_save / dt)
    tt = np.arange(0, Tmax, dt)  # array of times
    nts = len(tt)  # number of time-steps
    nts_save = nts - it_save  # number of time-steps to save

    x, y = octogon_positions()
    s_c = config['s_c']
    a, OmegaE = octogon_adjacency(s_c)
    I, J = np.where(np.triu(a))

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
        print('Creating GIF...')
        this_dir = os.path.dirname(os.path.abspath(__file__))
        input_directory = join(this_dir, 'animations')
        output_gif_name = "kuramoto_octagons.gif"
        create_gif(input_directory, output_gif_name)
    return R1a, R2a

import matplotlib.pyplot as plt

def plot_order_parameters(X1a, X2a, figsize=(8, 6), dpi=300):
    """
    Plots the magnitude and phase of two order parameters (R1a and R2a) over time.

    Args:
        X1a (numpy.ndarray): A 1D array containing the complex order parameter X1a values.
        X2a (numpy.ndarray): A 1D array containing the complex order parameter X2a values.
        figsize (tuple, optional): Size of the figure in inches. Defaults to (8, 6).
        dpi (int, optional): Resolution of the figure in dots per inch. Defaults to 300.
    """

    fig, ax = plt.subplots(2, 2, dpi=dpi, figsize=figsize)

    # Plot R1a
    ax[0, 0].plot(np.abs(X1a), linewidth=2)
    ax[0, 0].set_xlabel('$t$')
    ax[0, 0].set_ylabel('$R_1$')

    ax[0, 1].plot(np.angle(X1a), linewidth=2)
    ax[0, 1].set_xlabel('$t$')
    ax[0, 1].set_ylabel('Complex phase of  $X_1$')

    # Plot R2a
    ax[1, 0].plot(np.abs(X2a), linewidth=2)
    ax[1, 0].set_xlabel('$t$')
    ax[1, 0].set_ylabel('$R_2$')

    ax[1, 1].plot(np.angle(X2a), linewidth=2)
    ax[1, 1].set_xlabel('$t$')
    ax[1, 1].set_ylabel('Complex phase of  $X_2$')

    plt.tight_layout()


def plot_order_parameters_v2(X1a, X2a, figsize=(8, 6), dpi=300):
    """
    Plots the magnitude and phase of two order parameters (R1a and R2a) over time.

    Args:
        R1a (numpy.ndarray): A 1D array containing the order parameter R1a values.
        R2a (numpy.ndarray): A 1D array containing the order parameter R2a values.
        figsize (tuple, optional): Size of the figure in inches. Defaults to (8, 6).
        dpi (int, optional): Resolution of the figure in dots per inch. Defaults to 300.
    """

    fig, ax = plt.subplots(2, 2, dpi=dpi, figsize=figsize)

    # Plot R1a
    ax[0, 0].plot(np.abs(X1a), linewidth=2)
    ax[0, 0].set_xlabel('$t$')
    ax[0, 0].set_ylabel('$R_1$')

    ax[0, 1].plot(np.real(X1a),np.imag(X1a), linewidth=2)
    ax[0, 1].set_xlabel('$Re(X_1)$')
    ax[0, 1].set_ylabel(' $Im(X_1)$')

    # Plot R2a
    ax[1, 0].plot(np.abs(X2a), linewidth=2)
    ax[1, 0].set_xlabel('$t$')
    ax[1, 0].set_ylabel('$R_2$')

    ax[1, 1].plot(np.real(X2a),np.imag(X2a), linewidth=2)
    ax[1, 1].set_xlabel('$Re(X_2)$')
    ax[1, 1].set_ylabel(' $Im(X_2)$')

    plt.tight_layout()
