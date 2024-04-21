## Finite Element Project

This project consists of three distinct parts: pre-processing, processing, and post-processing of finite element problems. Each part was developed as part of the LEPL1110 course at UCLouvain on finite elements.

### Authors

- Mathis Delsart
- Adrien Antonutti

### Course Information

- **Course:** LEPL1110 - Finite Elements
- **Institution:** Université catholique de Louvain (UCLouvain)

### Project Structure

The project is organized into three main components:

1. **Pre-Processor**
2. **Finite Elements Processor**
3. **Finite Elements Post-Processor**

Each component is responsible for a specific step in the finite element analysis process.

### Build and Execution Instructions

To build and execute the entire project, you can use the provided automation script `main.c`. This script handles the construction and execution steps for each part of the project.
You can also build and execute each part separately by following the instructions in the respective directories.

#### Automation Script

The `main.c` script supports various options to customize the build and execution process. Here are some key options:

- `-m`: Disable mesh visualization.
- `-r`: Disable result visualization.
- `-e`: Enable the use of predefined examples.
- `-a`: Enable animation (only for post-processing).
- `-p`: Disable plotting (only for post-processing).
- `-h`: Display the help message.

### Directory Contents

- **Ressources**: Contains additional explanations on linear elasticity and instructions.
- **Rapport**: Contains the project report and its sources.

### Archive Validation

A script named `validate.py` is included to check the conformity of the archive. You can use it to ensure that all necessary files are present and properly organized.
This script come from the [\[Miguel-projects\](](https://github.com/MiguelDLC/femtools)