# Finite Elements Pre-Processor

## Authors

- **Student:** Mathis Delsart
- **Student:** Adrien Antonutti

## Course Info

- **Course:** LEPL1110
- **Institution:** UCLouvain

## Description

This repository contains the preprocessor developed for the LEPL1110 course at UCLouvain focusing on Finite Elements. The preprocessor serves the purpose of generating the mesh (saved in a `mesh.txt` file) and defining the problem (saved in a `problem.txt` file).

## Requirements

The preprocessor relies on the following libraries:

- [CMake](https://cmake.org/)
- [Gmsh](https://gmsh.info/)
- [OpenGl](https://www.opengl.org/)
- [GLFW](https://www.glfw.org/)

## Build Instructions

To build the preprocessor, execute the following commands in your terminal:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

To run the preprocessor, use the following command:

```bash
./myFem [-m] [-s] [-u] [-b] [-h]
```

### Options:

- `-m`: Disable the mesh visualizer.
- `-s`: Start the program with the bridge without stay cables and pylon.
- `-u`: Start the program with the U example mesh.
- `-b`: Start the program with the beam mesh.
- `-h`: Display this help message.

## Mesh Types

You can choose between triangle and quadrilateral mesh types by modifying the `elementType` variable in the `main.c` file:

- `FEM_TRIANGLE`: Triangle mesh.
- `FEM_QUAD`: Quadrilateral mesh.

## Elastic Cases

Choose between different types of problems by modifying the theCase variable in the main.c file:

- `PLANAR_STRESS`: Plane stress problem.
- `PLANAR_STRAIN`: Plane strain problem.
- `AXISYM`: Axisymmetric problem.


## Example Paths

By default, the program creates the complete bridge geometry with the PLANAR_STRESS case for the problem.