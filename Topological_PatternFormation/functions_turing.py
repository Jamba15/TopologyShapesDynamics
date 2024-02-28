import os
import sys
file_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(file_dir)
sys.path.append(parent_dir)
import fnmatch
import os
from os.path import join, dirname, abspath
import yaml
import numpy as np

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


class PrettySafeLoader(yaml.SafeLoader):
    def construct_python_tuple(self, node):
        return tuple(self.construct_sequence(node))

    # create constructor for !path tag
    def construct_path(self, node):
        return os.path.normpath(self.construct_scalar(node))



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


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):

        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


""" Function for dynamical system integration"""
def initial_condition(nodes, links, config, pert=1E-4):
    """
    This function defines the initial conditions for the dynamical system
    :param nodes: Dimension of the u vector (number of nodes)
    :param links:  Dimension of the v vector (number of links)
    :param config: configuration dictionary containing the "fixed_points" subdictionary
    :param pert: Perturbation strength
    :return: The vector of initial conditions concatenated as (u0, v0)
    """

    u_star = config['fixed_points']['u_star']
    v_star = config['fixed_points']['v_star']
    # Define the initial conditions with random perturbations
    u_pert = 1 - 2 * np.random.rand(nodes, 1)  # Random array for u perturbation
    v_pert = 1 - 2 * np.random.rand(links, 1)  # Random array for v perturbation

    # Assuming u_star and v_star are the fixed points for u and v, respectively
    u0 = u_star + pert * u_pert.flatten()  # Flatten to convert from 2D to 1D if necessary
    v0 = v_star + pert * v_pert.flatten()  # Flatten to convert from 2D to 1D if necessary

    # Combine u0 and v0 to create the initial condition array for the entire system
    y0 = np.concatenate((u0, v0))

    t_ini = config['integration_parameters']['t_ini']
    t_final = config['integration_parameters']['t_final']
    return t_ini, t_final, y0


# Define the system of differential equations to integrate
def system_to_integrate(t, y, config, B, L0, L1):
    # Define and parse the model_parameters
    a = config['model_parameters']['a']
    alpha = config['model_parameters']['alpha']
    b = config['model_parameters']['b']
    beta = config['model_parameters']['beta']
    c = config['model_parameters']['c']
    gamma = config['model_parameters']['gamma']

    # Define and parse the diffusion_coefficients
    D0 = config['diffusion_coefficients'].get('D0', 0)
    D1 = config['diffusion_coefficients'].get('D1', 0)
    D01 = config['diffusion_coefficients'].get('D01', 0)
    D10 = config['diffusion_coefficients'].get('D10', 0)

    # The vector y contains both u and v variables
    # Let's split y into u and v. Assuming u and v are concatenated one after the other in y:
    num_nodes = np.shape(L0)[0]  # Assuming L0 is a square matrix and num_nodes is the number of nodes

    u = y[:num_nodes]
    v = y[num_nodes:]

    # Compute f and g as defined in the MATLAB code
    f = -a * u - b * u**3 + c * (B @ v)
    g = -alpha * v - beta * v**3 + gamma * (B.T @ u)

    # Compute the derivatives of u and v
    du = f - D0 * (L0 @ u) - D01 * L0 @ (B @ v)
    dv = g - D1 * (L1 @ v) - D10 * L1 @ (B.T @ u)

    # Return the concatenated derivatives of u and v
    return np.concatenate((du, dv))

