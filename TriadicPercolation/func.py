import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix, find
import matplotlib.pyplot as plt


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
  a = a < (c/(N-1))
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

  for p in np.arange(0, 1, 0.02):

    R = np.zeros([Tmax, 1])   # Store the order parameter R for each time step.
    print(p)
    pL0 = 0.11    # Initial condition
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
        s[node] = 1   # Only nodes in the giant component are active.

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
    R_list = np.append(R_list, R[len(R)-num:len(R)])
  return p_list, R_list


def triadic_percolation_theory_poisson(c, cp, cn, Tmax, num):
    '''
      Theoretical calculation of the orbit diagram of the dynamic

      c: Average structural degree of the Poisson network
      cp: Average positive regulatory degree of the Poisson regulatory network
      cn: Average negative regulatory degree of the Poisson regulatory network
      Tmax: max duration of the dynamic
      num: Number of sample used to plot the orbit diagram
    '''
    R_list = []
    p_list = []

    for p in np.arange(0, 1, 0.01):
        R = []
        SN2 = 0.9
        pL = 0.14
        for i in range(Tmax):
            for nrun in range(100):
                SN2 = (1 - np.exp(-c * pL * SN2))
            pL = p * np.exp(-cn * SN2) * (1 - np.exp(-cp * SN2))
            R.append(SN2)

        p_list = np.append(p_list, p * np.ones([num, 1]))
        R_list = np.append(R_list, R[len(R)-num:len(R)])
    return p_list, R_list



