#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

int main(int argc, char *argv[])
{
    bool meshVisualizer = true;
    bool resultVisualizer = true;
    bool exampleUsage = false;
    bool animation = false;
    
    int opt;
    while ((opt = getopt(argc, argv, "aemrvh")) != -1)
    {
        switch (opt)
        {
            case 'a':
                animation = true;
                break;
            case 'e':
                exampleUsage = true;
                break;
            case 'm':
                meshVisualizer = false;
                break;
            case 'r':
                resultVisualizer = false;
                break;
            case 'v':
                meshVisualizer = false;
                resultVisualizer = false;
                break;
            case 'h':
                printf("Usage: %s [-e] [-v] [-r] [-a]\n", argv[0]);
                printf("Options:\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -m : Disable the mesh visualizer\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -v : Disable both visualizers\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-e] [-v] [-r] [-a] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }
    
    char *nameDirectory[] = {"Project/PreProcessing/build/", "../../Processing/build/", "../../PostProcessing/build/"};
    char *stages[] = {"Pre-Processing", "Processing", "Post-Processing"};
    int nbProgram = 3;

    for (int i = 0; i < nbProgram; i++)
    {
        if (access(nameDirectory[i], F_OK) == -1)
        {
            // Create the directory 'build' if it does not exist
            char command[50];
            sprintf(command, "mkdir -p %s", nameDirectory[i]);
            system(command);
        }
        
        // Change the current directory to the 'build' directory
        if (chdir(nameDirectory[i]) != 0) { perror("Erreur lors du changement de répertoire"); return EXIT_FAILURE; }

        // Execute the command "cmake .." to generate the Makefile
        system("cmake .. > /dev/null 2>&1");

        // Print the current step
        printf("\n\n/**************************/\n");
        printf("Stage %d : %s\n", i + 1, stages[i]);
        printf("/**************************/\n\n");

        // Execute the command "make" to build the program "./myFem"
        int ret_make = system("make 2>&1");
        if (ret_make != 0)
        {
            fprintf(stderr, "Erreur lors de l'exécution de la commande make lors de l'étape '%s'\n", stages[i]);
            exit(EXIT_FAILURE);
        }

        // Define the command to execute
        char command[40] = "./myFem";
        if (i == 0 && !meshVisualizer)   { strcat(command, " -m"); }
        if (i == 2 && !resultVisualizer) { strcat(command, " -r"); }
        if (i != 0 && animation)         { strcat(command, " -a"); }
        if (exampleUsage)                { strcat(command, " -e"); }
        strcat(command, " 2>&1");

        printf("\n\nStep %d : %s\n", i + 1, stages[i]);
        printf("Commande : %s\n\n", command);

        int ret_myFem = system(command);
        if (ret_myFem != 0)
        {
            fprintf(stderr, "Erreur lors de l'exécution de la commande '%s' lors de l'étape '%s'\n", command, stages[i]);
            exit(EXIT_FAILURE);
        }
    }

    return EXIT_SUCCESS;
}