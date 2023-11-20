import os
from scipy.io import loadmat
from numpy import savez
from os.path import dirname, abspath

# Estrai la directory del file corrente
file_directory = dirname(abspath(__file__))
os.chdir(file_directory)

if not os.path.exists('HyperCube_Boundary.m'):
    raise FileNotFoundError('Torus_Boundary.m not found in current directory')

data = loadmat('HyperCube_Boundary.m')

# Adesso data è un dizionario con le chiavi corrispondenti ai nomi delle variabili salvate.
# Puoi accedere agli array NumPy come segue:
B1 = data['B1']
B2 = data['B2']
B3 = data['B3']
B4 = data['B4']

#%%
# Salva le matrici con 4x4 in un file .npz mentre le altre in un altro
savez('B_Hypercube.npz', B1=B1, B2=B2, B3=B3, B4=B4)
