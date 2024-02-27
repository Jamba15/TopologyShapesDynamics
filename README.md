
# TODO: 
- Explain better each subsection of Experiments section
- Change folders to something more related
- Abstract/Introductory section
- ...
  
# Topology Shapes Dynamics
Higher-order networks capture the many-body interactions present in various empirical complex systems, revealing novel phenomena that shed new light on the interplay between the systems' topology and dynamics. 
The emerging theory of higher-order topological dynamics, combining higher-order structures with discrete topology and non-linear dynamics, has the potential to play a fundamental role for the understanding of complex systems such as the brain and climate, as well as for the development of a new generation of AI algorithms.
A new theoretical framework to describe network dynamics that goes beyond the node-centered description adopted so far thus emerges. In this novel theoretical framework, the dynamics of a network is encoded by topological signals, i.e., variables associated to  the nodes as well as to 
the edges (like fluxes), to the triangles, or even to other higher-order cells of higher-order networks.
One important challenge is  to model, mine, and process these signals to formulate a deeper physical and mathematical theory of complex systems. Recently it has been shown that topological signals lead to the emergence of novel types of dynamical states and collective phenomena such as topological synchronization, topological  pattern formation and triadic percolation. They offer novel paradigms to understand  how topology shapes dynamics, how dynamics learns the underlying network topology, and how topology  varies dynamically. 
This Perspective aims at guiding physicists, mathematicians, computer scientists and network scientists from different disciplines, in the growing field of topological signals, as well as at delineating challenges that must be addressed by future research. 

## Repository Overview

This repository contains the code to reproduce the results presented in the paper "Topology Shapes Dynamics on Higher-Order Networks". The repository is structured as follows:

* Each folder contains the code to reproduce the results of a given paper cited in the main work.
* The papers are organized with:
    * A configuration file that contains the parameters necessary for the simulations.
    * A Python file that contains the functions.
    * A Jupyter Notebook that contains the code to be launched.
* The main directory of the repository also contains a `yaml,yml` file with the configurations to give to Anaconda to create the environment.
Inspect always the configuration file as it contatins also a description of the parameters used in the simulations.
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

 ### TopKuramoto:
 contains the code for Topological Kuramoto model for eddge topological signals over simplicial and cell complexes.
 
  ### TopGlobalSync:
 contains the code for Topological Global Synchronization over simplicial and 
 cell complex using  Stuart-Landau oscillators.

 ### TopTuringPattern:
 contains the code for Topological Turing patterns for node and edge topological signals.
 
 ### DiracSignalProcess:
 contains the code for Dirac Signal Processing of node, edge and triangle topological signals over simplicial complexes.

 ### TriadicPercolation:
 contains the code for Triadic percolation on Poisson networks with triadic interactions.


## Citation

If you use this code, please cite the following paper:

```
[Topology Shapes Dynamics of Higher Order Networks]
```
