## README

### Repository Overview

This repository contains the code to reproduce the results presented in the paper "Paper Title: [https://paper.co/](https://paper.co/)". The repository is structured as follows:

* Each folder contains the code to reproduce the results of a given paper cited in the main work.
* The papers are organized with:
    * A configuration file that contains the parameters necessary for the simulations.
    * A Python file that contains the functions.
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

7. Install the required packages by running the following command in the terminal:

```
pip install -r requirements.txt
```

### Usage

To reproduce the results of a given paper, follow these steps:

1. Open the Jupyter Notebook in the folder of the paper.
2. Change the parameters in the configuration file if necessary.
3. Run the cells of the Jupyter Notebook.

### Citation

If you use this code, please cite the following paper:

```
[Paper Title](https://paper.co/)
```

### Contact

If you have any questions, please contact [Your Name] at [Your Email].