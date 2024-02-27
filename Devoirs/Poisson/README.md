# Instructions to Compile this Homework (Poisson)

## Compilation Steps:

### 1. Create Gmsh Directory:
Make a directory named `gmsh` using the command `mkdir gmsh`. Download the Gmsh Sdk version at this link [Gmsh Software Development Kit](https://www.gmsh.info/) and place the Gmsh SDK version into this folder (don't forget to unzip the folder). The path should resemble `gmsh/.../lib` where "..." denotes the version of the Gmsh SDK.
For example, for MACOS ARM: `gmsh-4.12.2-MacOSARM-sdk`.
This directory should be placed in the root directory of the homework, ./Geomesh.

### 2. Create Build Directory:
Create a directory named `build` with the command `mkdir build`.

### 3. Run CMake:
Navigate into the `build` directory and run `cmake ..`.

### 4. Compile and Execute:
Run `make` inside the `build` directory to compile the project.
Finally, run `./myFem` to execute the compiled program.

## Error Handling:
If an error occurs during the compilation process, you can follow these steps to resolve it:

### 1. Delete Gmsh Directory:
If an error occurs, it's advisable to remove the `gmsh` directory and download the Gmsh SDK again. Use the following commands to remove the `gmsh` directory: `rm -rf gmsh`

### 2. Delete Build Directory:
If an error occurs, it's often helpful to start fresh by deleting the build directory. Use the following command to remove the build directory: `rm -rf build`

### 3. Retry Compilation:
Once the `build` and `gmsh` directories have been removed, you can attempt the compilation process again by following the previous instructions.

By removing the `build` and `gmsh` directories, you can eliminate any potential conflicts or issues that may have arisen during the compilation process.