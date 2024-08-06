# Topology Shapes Dynamics
Higher-order networks encode complex interactions within systems, offering deeper insights into the relationship between
structure and behavior. By integrating higher-order structures, topology, and dynamics, we can develop powerful new models 
for systems ranging from neuroscience to climate science. Network dynamics can be described using topological signals – variables 
associated with nodes, edges, and higher-order network elements. Understanding these signals will reveal new dynamical states
(e.g., topological synchronization) and advance our understanding of how networks evolve. This repository contains the code to
reproduce the results presented in the paper "Topology Shapes Dynamics on Higher-Order Networks" and shows numerically how
topological signals shape the dynamics of higher-order networks. 

## Repository Overview

This repository contains the code to reproduce the results presented in the paper "Topology Shapes Dynamics on Higher-Order Networks". The repository is structured as follows:

* Each folder contains the code to reproduce the results of a given paper cited in the main work.
* The folders are organized with:
    * A configuration folder that contains the parameters necessary for the simulations. Always inspect the configuration file as it contatins also a description of the parameters used in the simulations.
    * Python files that contain functions or the main script.
    * A Jupyter Notebook that contains the code to be launched.
* The main directory of the repository also contains a `yaml` file with the configurations to give to Anaconda to create the environment.
### Installation

To install the environment, follow these steps:

1. Install Anaconda from [https://www.anaconda.com/](https://www.anaconda.com/).
2. Clone this repository to your local machine.
3. Open the `environment.yaml` file in a text editor.
4. Change the paths to the directories where you want to install the packages.
5. Create the environment by running the following command in the terminal:

```
conda env create -f environment.yaml
```

6. Activate the environment by running the following command in the terminal:

```
conda activate [environment_name]
```

### Usage

To reproduce the results of a given paper, follow these steps:

1. Open the Jupyter Notebook in the folder of the paper.
2. Change the parameters in the configuration file if necessary.
3. Run the cells of the Jupyter Notebook.

# Experiments 

 ### Topological_Kuramoto:
 Contains the code for **Topological Kuramoto model** for edge topological signals over simplicial and cell complexes. The 
code refers to Supplementary Material of the paper, *Section II D1*.
 
 ### Topological_GlobalSync:
 Contains the code for **Topological Global Synchronization** over simplicial and cell complex using Stuart-Landau oscillators.
The code refers to Supplementary Material of the paper, *Section II D2*.

 ### Topological_PatternFormation:
 Contains the code for **Topological Turing patterns** for node and edge topological signals. The code refers to Supplementary
 Material of the paper, section III C.
 
 ### Topological_SignalProcess:
 Contains the code for **Dirac Signal Processing** of node, edge and triangle topological signals over simplicial complexes.
The code refers to Supplementary Material of the paper, *Section IV B*.
 ### TriadicPercolation:
 Contains the code for **Triadic Percolation** on Poisson networks with triadic interactions. The code refers to Supplementary
 Material of the paper, *Section V B*.

### Supplementary Videos
[![Supplementary Video 4]]([https://www.youtube.com/watch?v=YOUTUBE_VIDEO_ID_HERE](https://youtu.be/ovvEvuMcACQ
))

https://youtu.be/ovvEvuMcACQ

## Citation
If you use this code, please cite the following paper:

```
[Topology Shapes Dynamics of Higher Order Networks]
```

## Code Contributors
Ginestra Bianconi, Timoteo Carletti, Lorenzo Giambagli, Ana P. Millán, Riccardo Muolo, Hanlin Sun
