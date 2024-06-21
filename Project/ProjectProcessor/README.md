# Finite Elements Processor

## Authors

- **Student:** Mathis Delsart
- **Student:** Adrien Antonutti

## Course Info

- **Course:** LEPL1110
- **Institution:** UCLouvain

## Description

This folder contains the processor for the LEPL1110 course at UCLouvain on Finite Elements.
The processor is used for solving problems related to linear elasticity.

## Requirements

The processor requires the following library:

- [CMake](https://cmake.org/)

## Build Instructions

To build the processor, use the following command:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

To run the processor, use the following command:

```bash
./myFem [-s] [-u] [-b] [-t] [-a] [-x] [-h]
```

### Options:

- `-s`: Start the program with the bridge without stay cables and pylon.
- `-u`: Start the program with the U example mesh.
- `-b`: Start the program with the beam mesh.
- `-t`: Time the program execution.
- `-a`: Launch the program 50 times to create an animation of the bridge's deformation (in ProjectPostProcessor).
- `-x`: Launch the program 50 times to create an animation of the car's movement (in ProjectPostProcessor).
- `-h`: Display this help message.

## Solver Types

The program supports different solver types, which can be selected by modifying the `typeSolver` variable in the `main.c` file:

- `FEM_FULL`: Full solver.
- `FEM_BAND`: Banded solver.

## Mesh Types

You can choose between triangle and quadrilateral mesh types by modifying the `elementType` variable in the `main.c` file:

- `FEM_TRIANGLE`: Triangle mesh.
- `FEM_QUAD`: Quadrilateral mesh.

## Renumbering Types

The renumbering type can be selected by modifying the `renumType` variable in the `main.c` file:

- `FEM_NO`: No renumbering.
- `FEM_XNUM`: Renumbering based on X coordinate.
- `FEM_YNUM`: Renumbering based on Y coordinate.
- `FEM_RCMK`: Reverse Cuthill-McKee renumbering.

## Example Paths

By default, the program reads mesh and problem files from the following paths:

- Mesh: `../data/mesh.txt`
- Problem: `../data/problem.txt`
- Solution: `../data/UV.txt`

You can modify these paths based on your specific requirements in the `main` function of `main.c`.