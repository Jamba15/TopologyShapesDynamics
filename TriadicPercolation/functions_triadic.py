import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix, find
import matplotlib.pyplot as plt
from tqdm import tqdm
from celluloid import Camera
from os.path import join, dirname, abspath
import yaml
from os import makedirs

current_folder = dirname(abspath(__file__))


def config_parser(config_file):
    # Look for the configuration file in the folder configurations
    config_file = join(current_folder, 'configurations', config_file)
    with open(config_file, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return config


def triadic_percolation_theory_poisson(config: dict, regulation="pos_neg"):
    '''
      Theoretical calculation of the orbit diagram of the dynamic

      c: Average structural degree of the Poisson network
      cp: Average positive regulatory degree of the Poisson regulatory network
      cn: Average negative regulatory degree of the Poisson regulatory network
      Tmax: max duration of the dynamic
      num: Number of sample used to plot the orbit diagram
      regulation: 
              "pos_neg": both positive and negative regulations are present
              "exclu_neg": Only negative regulations are present and positive regulations are not required
              "none": simple percolation without regulatory interactions
    '''

    c, cp, cn, Tmax, num = config['c'], config['cp'], config['cn'], config['Tmax'], config['num']

    R_list = []
    p_list = []

    for p in tqdm(np.arange(0, 1, 0.01)):
        R = []
        SN2 = 0.9
        pL = 0.04
        for i in range(Tmax):
            for nrun in range(100):
                SN2 = (1 - np.exp(-c * pL * SN2))
            if regulation == "pos_neg":
                pL = p * np.exp(-cn * SN2) * (1 - np.exp(-cp * SN2))
            elif regulation == "exclu_neg":
                pL = p * np.exp(-cn * SN2)
            elif regulation == "none":
                pL = p
            R.append(SN2)

        p_list = np.append(p_list, p * np.ones([num, 1]))
        R_list = np.append(R_list, R[len(R) - num:len(R)])
    return p_list, R_list


def triadic_percolation_simulation_poisson_movie(config, movie_name, regulation="pos_neg"):
    '''
    Monte Carlo simulation of the orbit diagram of the dynamic

    c: Average structural degree of the Poisson network
    cp: Average positive regulatory degree of the Poisson regulatory network
    cn: Average negative regulatory degree of the Poisson regulatory network
    Tmax: max duration of the dynamic
    num: Number of sample used to plot the orbit diagram
    regulation: 
              "pos_neg": both positive and negative regulations are present
              "exclu_neg": Only negative regulations are present and positive regulations are not required

  '''
    N, p, c, cp, cn, Tmax = config['N'], config['p'], config['c'], config['cp'], config['cn'], config['Tmax']
    a = np.random.rand(N, N)
    a = a < (c / (N - 1))
    a = a.astype(int)
    a = np.triu(a, 1)
    a = a + a.transpose()
    [I, J, V] = find(np.triu(a))  # Find all the links of the adjacency matrix

    '''
    adj_pos[i,l]=1 if node i positively regulates link l
    adj_pos[i,l]=1 if node i negatively regulates link l
  '''

    L = sum(sum(a)) / 2  # Number of links
    rand_reg = np.random.rand(N, L.astype(int))
    adj_pos = (rand_reg < cp / N)
    adj_neg = (rand_reg > cp / N) * (rand_reg < (cp + cn) / N)

    R = np.zeros([Tmax + 1, 1])  # Store the order parameter R for each time step.
    pL0 = 0.04  # Initial condition
    xL = np.random.rand(L.astype(int), 1)

    # Initialize links. state[l]=1 if the link l is active.
    state = (xL < pL0)
    state = state.astype(float)

    [linkID, _, _] = find(state)  # return the list of active links
    new_V = [V[i] for i in linkID]
    new_I = [I[i] for i in linkID]
    new_J = [J[i] for i in linkID]

    # Construct a new network with only active links.
    adj = csr_matrix((new_V, (new_I, new_J)), shape=(N, N))
    adj = adj + np.transpose(adj)

    G = nx.from_numpy_array(adj)

    Gcc = max(nx.connected_components(G), key=len)  # Find the giantcomponent.
    G0 = G.subgraph(Gcc)

    # Initialize the state of nodes. s[i]=1 if node i is in the giant component.
    s = np.zeros([N, 1])
    if len(G0) > 1:
        for node in G0.nodes():
            s[node] = 1  # Only nodes in the giant component are active.

    R[0] = len(G0) / N
    pos = nx.random_layout(G, seed=1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    camera = Camera(fig)
    plt.ioff()
    for it in range(1, Tmax + 1):
        xL = np.random.rand(L.astype(int), 1)
        if regulation == "pos_neg":
            state = (xL < p) * ((adj_pos.transpose() @ s) > 0) * (
                        (adj_neg.transpose() @ s) == 0)  # Conditions for a link being active

        elif regulation == "exclu_neg":  # The case where positive regulations are not required, the regulation formed exclusively by negative regulations
            state = (xL < p) * ((adj_neg.transpose() @ s) == 0)  # Conditions for a link being active

        state = state.astype(float)
        [linkID, _, _] = find(state)
        new_V = [V[i] for i in linkID]
        new_I = [I[i] for i in linkID]
        new_J = [J[i] for i in linkID]
        adj = csr_matrix((new_V, (new_I, new_J)), shape=(N, N))
        adj = adj + np.transpose(adj)
        G = nx.from_scipy_sparse_array(adj)
        Gcc = max(nx.connected_components(G), key=len)
        G0 = nx.Graph()
        G0 = G.subgraph(Gcc)
        s = np.zeros([N, 1])
        if len(G0) > 1:
            for node in G0.nodes():
                s[node] = 1
        R[it] = len(G0) / N

        node_color_map = []
        edge_color_map = []

        outline_color = "black"
        active_node_color = (243 / 255, 175 / 255, 105 / 255, 1)
        edge_color = (243 / 255, 175 / 255, 105 / 255, 1)
        inactive_node_color = (251 / 255, 231 / 255, 211 / 255, 1)
        transparent = (1, 1, 1, 0)

        for node in G:
            if node in G0.nodes():
                node_color_map.append(active_node_color)
            else:
                node_color_map.append(inactive_node_color)

        for edge in G.edges():
            if edge in G0.edges():
                edge_color_map.append(edge_color)
            else:
                edge_color_map.append(transparent)

        plt.rcParams["font.family"] = "Times New Roman"
        nx.draw(G, pos, ax=ax1, node_size=200, edgecolors=outline_color, node_color=node_color_map,
                edge_color=edge_color_map, width=3)
        ax1.axis('off')
        ax2.plot(range(1, it + 1), R[1:it + 1], color=active_node_color, linewidth=3, marker='o')
        ax2.axis('on')
        ax2.set_xticks([0, 10, 20, 30, 40, 50])
        ax2.set_xticklabels(['$0$', '$10$', '$20$', '$30$', '$40$', '$50$'], fontsize=18)
        ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax2.set_yticklabels(['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'], fontsize=18)
        ax2.set_xlim(-0.2, 50.2)
        ax2.set_ylim(0, 1)
        ax2.set_xlabel("Time", fontsize=18)
        ax2.set_ylabel("Order parameter", fontsize=18)
        fig.tight_layout()
        camera.snap()

    animation = camera.animate()
    # Save animation in animation folder
    folder_name = join(current_folder, 'animations')
    # Create the folder if it does not exist
    makedirs(folder_name, exist_ok=True)

    animation.save(join(folder_name, movie_name+".gif"), writer='imagemagick', fps=5)
    return 0


def triadic_percolation_theory_poisson_time_series(config, regulation="pos_neg"):
    '''
      Theoretical calculation of the time series of the dynamic

      c: Average structural degree of the Poisson network
      p: probability of retaining a link after the regulation process
      cp: Average positive regulatory degree of the Poisson regulatory network
      cn: Average negative regulatory degree of the Poisson regulatory network
      Tmax: max duration of the dynamic
      num: Number of sample used to plot the orbit diagram
      regulation: 
              "pos_neg": both positive and negative regulations are present
              "exclu_neg": Only negative regulations are present and positive regulations are not required
              "none": simple percolation without regulatory interactions
    '''
    c, p, cp, cn, Tmax, num = config['c'], config['p'], config['cp'], config['cn'], config['Tmax'], config['num']
    R_list = []
    T_list = []
    SN2 = 0.9
    pL = 0.04
    for i in range(Tmax):
        for nrun in range(100):
            SN2 = (1 - np.exp(-c * pL * SN2))
        if regulation == "pos_neg":
            pL = p * np.exp(-cn * SN2) * (1 - np.exp(-cp * SN2))
        elif regulation == "exclu_neg":
            pL = p * np.exp(-cn * SN2)
        elif regulation == "none":
            pL = p
        if i > Tmax - num:
            R_list.append(SN2)
            T_list.append(i)

    return T_list, R_list


def plot_orbit_diagram(p_list, R_list):
    plt.clf()
    plt.figure()
    active_node_color = (243 / 255, 175 / 255, 105 / 255, 1)

    plt.scatter(p_list, R_list, color=active_node_color, linewidths=0.5, marker='.')
    plt.xticks(ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'],
               fontsize=18)
    plt.yticks(ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'],
               fontsize=18)

    plt.xlabel("$p$", fontsize=18)
    plt.ylabel("$R$", fontsize=18)
    plt.show()


def plot_time_series(T_list, R_list):
    plt.clf()
    plt.figure()
    active_node_color = (243 / 255, 175 / 255, 105 / 255, 1)
    plt.plot(T_list, R_list, color=active_node_color)
    plt.yticks(ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'],
               fontsize=18)
    plt.xticks(fontsize=18)
    plt.xlabel("$T$", fontsize=18)
    plt.ylabel("$R$", fontsize=18)
    plt.show()


def triadic_percolation_simulation_poisson(N, c, cp, cn, Tmax, num):
    '''
    Monte Carlo simulation of the orbit diagram of the dynamic

    c: Average structural degree of the Poisson network
    cp: Average positive regulatory degree of the Poisson regulatory network
    cn: Average negative regulatory degree of the Poisson regulatory network
    Tmax: max duration of the dynamic
    num: Number of sample used to plot the orbit diagram
  '''

    a = np.random.rand(N, N)
    a = a < (c / (N - 1))
    a = a.astype(int)
    a = np.triu(a, 1)
    a = a + a.transpose()
    [I, J, V] = find(np.triu(a))  # Find all the links of the adjacency matrix

    '''
    adj_pos[i,l]=1 if node i positively regulates link l
    adj_pos[i,l]=1 if node i negatively regulates link l
  '''

    L = sum(sum(a)) / 2  # Number of links
    rand_reg = np.random.rand(N, L.astype(int))
    adj_pos = (rand_reg < cp / N)
    adj_neg = (rand_reg > cp / N) * (rand_reg < (cp + cn) / N)

    p_list = []
    R_list = []

    for p in tqdm(np.arange(0, 1, 0.02)):

        R = np.zeros([Tmax, 1])  # Store the order parameter R for each time step.
        pL0 = 0.04  # Initial condition
        xL = np.random.rand(L.astype(int), 1)
        # Initialize links. state[l]=1 if the link l is active.
        state = (xL < pL0)
        state = state.astype(float)

        [linkID, _, _] = find(state)  # return the list of active links
        new_V = [V[i] for i in linkID]
        new_I = [I[i] for i in linkID]
        new_J = [J[i] for i in linkID]

        # Construct a new network with only active links.
        adj = csr_matrix((new_V, (new_I, new_J)), shape=(N, N))
        adj = adj + np.transpose(adj)
        G = nx.from_scipy_sparse_array(adj)
        Gcc = max(nx.connected_components(G), key=len)  # Find the giantcomponent.
        G0 = G.subgraph(Gcc)

        # Initialize the state of nodes. s[i]=1 if node i is in the giant component.
        s = np.zeros([N, 1])
        if len(G0) > 1:
            for node in G0.nodes():
                s[node] = 1  # Only nodes in the giant component are active.

        R[0] = len(G0) / N

        for it in range(1, Tmax):
            xL = np.random.rand(L.astype(int), 1)
            state = (xL < p) * ((adj_pos.transpose() @ s) > 0) * \
                    ((adj_neg.transpose() @ s) == 0)  # Conditions for a link being active
            state = state.astype(float)
            [linkID, _, _] = find(state)
            new_V = [V[i] for i in linkID]
            new_I = [I[i] for i in linkID]
            new_J = [J[i] for i in linkID]
            adj = csr_matrix((new_V, (new_I, new_J)), shape=(N, N))
            adj = adj + np.transpose(adj)
            G = nx.from_scipy_sparse_array(adj)
            Gcc = max(nx.connected_components(G), key=len)
            G0 = G.subgraph(Gcc)

            s = np.zeros([N, 1])
            node_list = np.array(G0.nodes())
            if len(G0) > 1:
                s[node_list.tolist()] = 1

            R[it] = len(G0) / N
        p_list = np.append(p_list, p * np.ones([num, 1]))
        R_list = np.append(R_list, R[len(R) - num:len(R)])
    return p_list, R_list
