import fnmatch
import os
from os.path import join, dirname, abspath
import yaml
import numpy as np
import re


def import_matrices(file_name):
    """
    This function imports the matrices from the .npz file
    :param file_name: The name of the file containing the matrices
    :return: The adjacency and boundary matrices
    """
    file_name = find(file_name, join(dirname(abspath(__file__)), 'Matrices'))[0]
    matrices = np.load(file_name)
    # Crea un nuovo dizionario per gli array modificati
    modified_matrices = {}

    # Itera attraverso gli array, convertendoli e salvandoli nel nuovo dizionario
    for key in matrices:
        modified_matrices[key] = matrices[key].astype(np.float32)

    return modified_matrices


def get_Laplacian(k, f_name='B_Torus.npz'):

    Bik = import_matrices(f_name)
    # Adjust for 0-indexing in Python

    # Handle Bk
    if k == 0:
        Bk = np.zeros_like(Bik['B1'])  # Assuming shape from Bik['B1']
    else:
        Bk = Bik[f'B{k}']

    # Handle Bkp1
    if k == 3:
        Bkp1 = np.zeros_like(Bik['B1'])  # Assuming shape from Bik['B1']
    else:
        Bkp1 = Bik[f'B{k + 1}']

    Lk = -(np.dot(Bk.T, Bk) + np.dot(Bkp1, Bkp1.T))

    Lkup = -np.dot(Bkp1, Bkp1.T)
    Lkdwn = -np.dot(Bk.T, Bk)

    return Lk, Lkup, Lkdwn


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


""" Function for dynamical system integration"""


def initial_condition(k, config):
    """
    This function defines the initial conditions for the dynamical system
    :param nodes: Dimension of the u vector (number of nodes)
    :param links:  Dimension of the v vector (number of links)
    :param config: configuration dictionary containing the "fixed_points" subdictionary
    :param pert: Perturbation strength
    :return: The vector of initial conditions concatenated as (u0, v0)
    """
    Lk, _, _ = get_Laplacian(k)
    nk = Lk.shape[0]

    # Define the initial conditions for the dynamical system
    w_init = np.random.rand(nk) + 1j * np.random.rand(nk)

    t_ini = config['integration_parameters']['t_ini']
    t_final = config['integration_parameters']['t_final']
    return t_ini, t_final, w_init, Lk


# Define the system of differential equations to integrate
def system_to_integrate(t, w, config, Lk):
    # Define and parse the model_parameters

    sigma = config['model_parameters']['sigma']
    beta = config['model_parameters']['beta']
    mu = config['model_parameters']['mu']
    m = config['model_parameters']['m']

    f = sigma * w - beta * np.abs(w) ** 2 * w
    Hz = w * np.abs(w) ** (m - 1)

    dw = f + mu * np.dot(Lk, Hz)

    # Return the derivative of w
    return dw


