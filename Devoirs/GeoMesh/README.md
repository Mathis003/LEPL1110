### Instructions to Compile this Homework (Geomesh)

## To compile this homework, follow these steps:

# 1. Create Gmsh Directory:
Make a directory named `gmsh` using the command `mkdir gmsh`. Download the Gmsh Sdk version at this link [Gmsh Software Development Kit](https://www.gmsh.info/) and place the Gmsh SDK version into this folder (don't forget to unzip the folder). The path should resemble `gmsh/.../lib` where "..." denotes the version of the Gmsh SDK.
For example, for MACOS ARM: `gmsh-4.12.2-MacOSARM-sdk`.
This directory should be placed in the root directory of the homework, ./Geomesh.

# 2. Create Build Directory:
Create a directory named `build` with the command `mkdir build`.

# 3. Run CMake:
Navigate into the `build` directory and run `cmake ..`.

# 4. Compile and Execute:
Run `make` inside the `build` directory to compile the project.
Finally, run `./myFem` to execute the compiled program.