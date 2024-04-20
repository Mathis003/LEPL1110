# Finite Elements Post-Processor

## Authors

- **Mathis Delsart**
- **Adrien Antonutti**

## Course Information

- **Course:** LEPL1110 - Finite Elements
- **Institution:** Université catholique de Louvain (UCLouvain)

## Description

This repository contains the postprocessor developed for the LEPL1110 course at UCLouvain focusing on Finite Elements. The postprocessor is designed for visualizing and analyzing the results of finite element analysis for structures like bridges.

## Requirements

The postprocessor requires the following libraries:

- [CMake](https://cmake.org/)
- [OpenGL](https://www.opengl.org/)
- [GLFW](https://www.glfw.org/)

## Build Instructions

To build the postprocessor, follow these steps:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

To run the postprocessor, use the following command:

```bash
./myFem [-r] [-p] [-a] [-u] [-b] [-s] [-h]
```

### Options:

- `-r`: Disable the result visualizer.
- `-p`: Disable the plot.
- `-a`: Enable the animation.
- `-u`: Start the program with the example UForm.
- `-b`: Start the program with the example beam.
- `-s`: Start the program with the example simplified.
- `-h`: Display this help message.

## Solver Types

The program supports different solver types, which can be selected by modifying the `typeSolver` variable in the `main.c` file:

- `FEM_FULL`: Full solver.
- `FEM_BAND`: Banded solver.

## Renumbering Types

The renumbering type can be selected by modifying the `renumType` variable in the `main.c` file:

- `FEM_NO`: No renumbering.
- `FEM_XNUM`: Renumbering based on X coordinate.
- `FEM_YNUM`: Renumbering based on Y coordinate.
- `FEM_RCMK`: Reverse Cuthill-McKee renumbering.

## Program Options

- `V`: Mesh and displacement norm.
- `D`: Domains.
- `N`: Next domain highlighted.
- `S`: Display the matrix.
- `X`: Display the horizontal forces.
- `Y`: Display the vertical forces.
- `U`: Display the sigmaXX (Only for complete bridge structure).
- `I`: Display the sigmaYY (Only for complete bridge structure).
- `O`: Display the sigmaXY (Only for complete bridge structure).
- `P`: Display the von Mises (Only for complete bridge structure).

## Additionnal scripts

# beam_plot.py

Plot the results of the beam example with the analytical solution plotted in the middle of the beam.

# constrains_plot.py

Plot the results of the constraints on the complete bridge structure.

# convergence_plot.py

Plot the convergence of the FEM method with the beam example.

# solver_plot.py

Plot the results of the solver on the complete bridge structure.

# plot.py

Plot the results of a structure with various parameters. Check the help message for more information.