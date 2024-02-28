import os
from scipy.io import loadmat
from numpy import savez
from os.path import dirname, abspath

# Estrai la directory del file corrente
file_directory = dirname(abspath(__file__))
os.chdir(file_directory)

if not os.path.exists('myNetworks.mat'):
    raise FileNotFoundError('myNetworks.mat not found in current directory')

data = loadmat('myNetworks.mat')

# Adesso data è un dizionario con le chiavi corrispondenti ai nomi delle variabili salvate.
# Puoi accedere agli array NumPy come segue:
A_2D_4x4 = data['A_2D_4x4']
B_2D_4x4 = data['B_2D_4x4']
A_benchmark = data['A_benchmark']
B_benchmark = data['B_benchmark']
#%%
# Salva le matrici con 4x4 in un file .npz mentre le altre in un altro
savez('Matrices_4x4.npz', adjacency=A_2D_4x4, boundary=B_2D_4x4)
savez('Matrices_benchmark.npz', adjacency=A_benchmark, boundary=B_benchmark)
