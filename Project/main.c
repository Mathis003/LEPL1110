#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

int main(int argc, char *argv[])
{
    bool meshVisualizer = true;
    bool resultVisualizer = true;

    // Deal with the options arguments
    int opt;
    while ((opt = getopt(argc, argv, "mrah")) != -1)
    {
        switch (opt)
        {
            case 'm':
                meshVisualizer = false;
                break;
            case 'r':
                resultVisualizer = false;
                break;
            case 'a':
                meshVisualizer = false;
                resultVisualizer = false;
                break;
            case 'h':
                printf("Usage: %s [-v] [-r] [-a]\n", argv[0]);
                printf("Options:\n");
                printf("  -m : Disable the mesh visualizer\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -a : Disable both visualizers\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-v] [-r] [-a] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }
    
    FILE *fp;
    const int MAXLINE = 999999;
    char res[MAXLINE];

    char *nameDirectory[] = {"Project/PreProcessing/build", "../../Processing/build", "../../PostProcessing/build"};
    char *stages[] = {"Pre-Processing", "Processing", "Post-Processing"};
    int nbProgram = 3;

    for (int i = 0; i < nbProgram; i++)
    {   
        // Change the directory to the current program
        if (chdir(nameDirectory[i]) != 0) { perror("Erreur lors du changement de répertoire"); return EXIT_FAILURE; }
        
        // Print the current step
        printf("\n\n/**************************/\n");
        printf("Stage %d : %s\n", i + 1, stages[i]);
        printf("/**************************/\n\n");

        // Define the command to execute
        char command[100] = "./myFem";
        if      (i == 0 && !meshVisualizer)   { sprintf(command, "./myFem -m"); }
        else if (i == 2 && !resultVisualizer) { sprintf(command, "./myFem -r"); }

        // Execute the command "make" to build the program "./myFem"
        system("make > /dev/null");

        // Execute the program "./myFem" with the option(s)
        fp = popen(command, "r");
        if (fp == NULL) { printf("Error : fp == NULL\n"); return EXIT_FAILURE; }

        // Display the output of the command
        while (fgets(res, MAXLINE, fp) != NULL) { printf("%s", res); }

        // Close the file
        if (pclose(fp) == -1) { printf("Error : pclose(fp) == -1\n"); return EXIT_FAILURE; }
    }

    return EXIT_SUCCESS;
}