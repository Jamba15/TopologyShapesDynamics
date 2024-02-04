import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors
from pylab import rcParams
import seaborn as sns
from os.path import join
sns.set_style("whitegrid", {'axes.grid': False})
import numpy as np
import networkx as nx
from scipy.sparse import csc_matrix
import os
import yaml

current_folder = os.path.dirname(os.path.abspath(__file__))

def config_parser(config_file):
    # Look for the configuration file in the folder configurations
    config_file = join(current_folder, 'configurations', config_file)
    with open(config_file, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return config



def process(D, s_input, m, gamma):
    Q_op = gamma * (D.dot(D) - 2 * m * D + m ** 2 * np.identity(np.shape(D)[0]))
    op = np.identity(np.shape(D)[0]) + Q_op
    return np.linalg.inv(op).dot((s_input))


def process_iterative(D, s_input, m, gamma, Nrun_max):
    Q_op = gamma * (D.dot(D) - 2 * m * D + m ** 2 * np.identity(np.shape(D)[0]))
    s_estimated = s_input
    for _ in np.linspace(1, Nrun_max, Nrun_max):
        s_estimated = Q_op.dot(s_estimated)
    return s_estimated


def optimize_m(D, s_noisy, s_true, m, gamma, tol):
    epsilon = 0.3
    list_m = []
    list_error_m = []
    list_it = []
    it = 0
    m_old = 10000
    s_est_m = None

    while (np.abs(m_old - m)) > tol:
        m_old = m

        it = it + 1
        s_est_m = process(D, s_noisy, m, gamma)
        m = (1 - epsilon) * m_old + epsilon * ((s_est_m.transpose()).dot(D)).dot(s_est_m) / (
            s_est_m.transpose().dot(s_est_m))
        error_m = np.linalg.norm(s_est_m - s_true, 2)

        list_m.append(m)
        list_error_m.append(error_m)
        list_it.append(it)

    return list_m, list_error_m, list_it, m, s_est_m


class signal_elaboration:
    def __init__(self, D):
        self.D = D
        w, v = np.linalg.eig(D)
        w_real = np.real_if_close(w)
        w_min = w_real[(w_real > 0) & (~np.isclose(w_real, 0))].min()
        self.m_min_index = np.where(w_real == w_min)[0]
        w_max = w_real[(w_real > 0) & (~np.isclose(w_real, 0))].max()
        self.m_max_index = np.where(w_real == w_max)[0]
        mask_neg_and_nonzero = (w_real < 0) & (~np.isclose(w_real, 0))
        self.v_neg = v[:, mask_neg_and_nonzero]
        self.w_real = w_real
        self.v = v

    def signal_processing(self, Gaussian, x_min, gamma, m_0, tol):
        alpha = 0.6
        s_1 = np.real(self.v[:, self.m_min_index[0]]) if x_min == True else np.real(self.v[:, self.m_max_index[0]])
        noise = np.random.normal(0, 1, s_1.shape[0]) if Gaussian == True else self.v_neg.dot(
            np.random.randn(self.v_neg.shape[1]))
        noise = (noise.real - s_1 * ((s_1.T).dot(noise.real)))
        noise = noise / np.linalg.norm(noise)

        s_noisy = s_1 + (1 - alpha) ** 0.5 * noise
        s_true = s_1

        list_m, list_error_m, list_it, m, s_est_m = optimize_m(self.D, s_noisy, s_true, m_0, gamma, tol)
        return list_m, list_error_m, list_it, m, s_est_m

    def visualization_plot(self, i1, Gaussian, x_min, lab, ax=None):
        # x_min=True selects true signal aligned to the eigenvector associated to the smallest positive eigenvalue
        # x_min= False selects true signal aligned to the eigenvector associated to the largest positive eigenvalue
        # Gaussian=True selects gaussian noise Gaussian=False select noise having only contribution coming from vector
        # with negative chiarality
        """
                Plots the results of signal processing with options for different noise types and eigenvalue alignments.

                Parameters:
                    i1 (int): Index for selecting the axis.
                    Gaussian (bool): Selects Gaussian noise if True, otherwise a different noise type.
                    x_min (bool): Aligns to the smallest/largest positive eigenvalue if True/False.
                    lab (bool): Determines if labels are included in the plot.
                    ax (matplotlib axis, optional): Axis to plot on. If None, a new figure is created.
                """

        if ax is None:
            fig, ax = plt.subplots(nrows=2, ncols=1)
            axis1, axis2 = ax[0], ax[1]
        else:
            axis1 = ax[0][i1]
            axis2 = ax[1][i1]

        gamma = 10.
        m_0 = 1.4 if x_min else 4.0
        tol = 0.001
        m_index = self.m_max_index[0] if not x_min else self.m_min_index[0]
        list_m, list_error_m, list_it, m, s_est_m = self.signal_processing(Gaussian, x_min, gamma, m_0, tol)

        true_m_line = self.w_real[m_index] * np.ones(np.shape(list_it)[0])
        if lab:
            axis1.plot(list_it, list_m, label='inferred m')
            axis1.plot(list_it, true_m_line, label='true m')
        else:
            axis1.plot(list_it, list_m)
            axis1.plot(list_it, true_m_line)

        axis2.plot(list_it, list_error_m, label='error')

        # Only call legend if labels are present
        if lab or len(axis2.get_legend_handles_labels()[0]) > 0:
            axis1.legend(frameon=False)
            axis2.legend(frameon=False)

        axis1.set_xlabel('$time$')
        axis1.set_ylabel('$m$')
        axis2.set_xlabel('$time$')
        axis2.set_ylabel('$\Delta s(m)/\Delta s(0)$')

    def plot_results(self, save_fig, **kwargs):
        rcParams['lines.linewidth'] = 2
        font = {'weight': 'bold',
                'size': 18}

        matplotlib.rc('font', **font)

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))

        self.visualization_plot(0, False, True, True, ax=axs)
        self.visualization_plot(1, False, False, False, ax=axs)

        plt.annotate('(a)', xy=(-1.15, 1.30), xycoords='axes fraction', size=18)
        plt.annotate('(b)', xy=(0.05, 1.30), xycoords='axes fraction', size=18)
        plt.annotate('(c)', xy=(-1.15, 0.05), xycoords='axes fraction', size=18)
        plt.annotate('(d)', xy=(0.05, 0.05), xycoords='axes fraction', size=18)
        plt.title('Signal Processing')
        plt.show()

        if save_fig:
            fig_name = kwargs.get('fig_name', 'Plot Results.pdf')
            fig.savefig(fig_name)


class SimplicialComplexConstructor:
    def __init__(self, **kwargs):
        data_folder = join(current_folder,kwargs.get('data_folder', 'Data'))
        csv_file = kwargs.get('csv_file', None)
        edge_list_file = kwargs.get('edge_list_file', None)
        triangle_file = kwargs.get('triangles_list_file', None)

        if csv_file is not None:
            df = pd.read_csv(join(data_folder, csv_file), sep=" ", header=None,
                             names=["layerID", "node1", "node2", "weight"])

            # Filter rows with layerID 1 and copy the relevant columns
            b1 = df[df['layerID'] == 1]
            edge_list1 = b1[['node1', 'node2']].copy()

            # Create a graph and add edges
            self.G = nx.Graph()
            self.G.add_edges_from(edge_list1.values)

            # Convert node labels to integers
            self.G = nx.convert_node_labels_to_integers(self.G, 1)
            self.edge_list = list(self.G.edges())

            # Calculate the number of nodes and length of the edge list
            self.N = self.G.number_of_nodes()
            self.L = len(self.edge_list)
            self.T = 0  # Number of triangles, set to 0 as per your original code

        if (triangle_file is not None) and (edge_list_file is not None):
            self.edge_list = pd.read_csv(join(data_folder, edge_list_file), sep=",", header=None, names=["node1", "node2"])

            # Read the triangle data
            dft = pd.read_csv(join(data_folder, triangle_file), sep=",", header=None, names=["node1", "node2", "node3"])
            self.triangles_list = np.array(list(dft.values))

            # Create a graph and add edges
            self.G = nx.Graph()
            self.G.add_edges_from(self.edge_list.values)

            # Calculating the number of nodes, edge list, and triangle list
            self.N = self.G.number_of_nodes()
            self.L = len(self.edge_list)
            self.edge_list = self.edge_list.values
            self.triangle_list = np.array([sorted(j) for j in self.triangles_list])
            self.T = len(self.triangle_list)

    def draw_graph(self):
        if not hasattr(self, 'triangle_list'):
        # Draw the graph with node labels
            pos = nx.spring_layout(self.G)  # Positioning of nodes
            nx.draw(self.G, pos, with_labels=True, node_color='lightblue', edge_color='gray')
            plt.show()

        else:
            plt.figure()
            plt.title("Simplex form edges and triangles", fontsize=20)
            pos = nx.spring_layout(self.G, iterations=100)
            nx.draw(self.G, pos)
            # fill the triangles
            for t in self.triangle_list:
                plt.gca().add_patch(mpatches.Polygon(np.array([pos[t[0]], pos[t[1]], pos[t[2]]]), alpha=0.2))

            plt.show()

    def create_B1_D1(self):
        # ASSUMES FIRST NODE IS INDEXED AS 1 AND EDGES ARE SORTED (n1,n2) with n1<n2
        num_edges = len(self.edge_list)
        data = [-1] * num_edges + [1] * num_edges
        row_ind = [e[0] - 1 for e in self.edge_list] + [e[1] - 1 for e in self.edge_list]
        col_ind = [i for i in range(len(self.edge_list))] * 2
        B1 = csc_matrix(
            (np.array(data), (np.array(row_ind), np.array(col_ind))), dtype=np.int8)

        B1 = B1.toarray()
        D1_top = np.concatenate([np.zeros((self.N, self.N)), B1], axis=1)
        D1_bottom = np.concatenate([B1.T, np.zeros((self.L, self.L))], axis=1)
        D1 = np.concatenate([D1_top, D1_bottom], axis=0)

        return B1, D1

    def create_B2_D2(self):
        if len(self.triangle_list) == 0:
            return csc_matrix([], shape=(len(self.edge_list), 0), dtype=np.int8)

        elist_dict = {tuple(sorted(j)): i for i, j in enumerate(self.edge_list)}

        data = []
        row_ind = []
        col_ind = []
        for i, t in enumerate(self.triangle_list):
            e1 = t[[0, 1]]
            e2 = t[[1, 2]]
            e3 = t[[0, 2]]

            data.append(1)
            row_ind.append(elist_dict[tuple(e1)])
            col_ind.append(i)

            data.append(1)
            row_ind.append(elist_dict[tuple(e2)])
            col_ind.append(i)

            data.append(-1)
            row_ind.append(elist_dict[tuple(e3)])
            col_ind.append(i)

        B2 = csc_matrix((np.array(data), (np.array(row_ind), np.array(
            col_ind))), shape=(len(self.edge_list), len(self.triangle_list)), dtype=np.int8)

        D2_top = np.concatenate([np.zeros((self.L, self.L)), B2.toarray()], axis=1)
        D2_bottom = np.concatenate([B2.toarray().T, np.zeros((self.T, self.T))], axis=1)
        D2 = np.concatenate([D2_top, D2_bottom], axis=0)
        return B2.toarray(), D2

    def create_gamma2(self): # TODO chiarire a cosa serve questa parte
        return np.block([[-np.eye(self.L),np.zeros((self.L,self.T))],[np.zeros((self.T,self.L)),np.eye(self.T)]])

