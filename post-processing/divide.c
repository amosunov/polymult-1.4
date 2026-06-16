/*
Takes as input a prefix of input binary files, divides their contents (unsigned integers) by a fixed positive unsigned integer, and saves the resulting values into an output file

EXAMPLE: ./sum 32 ~/res ~/a 3

The input files have names a0, a1, ..., a31
The output files are called res0, res1, ..., res31
If we view both files as huge vectors, then res0=(1/3)*a0, res1(1/3)*a1, ..., res31=(1/3)*a31 

Notice that all files must store binary unsigned int values, have exactly the same nonzero size, and this size must be divisible by BUFSIZE

Compile with clang -fopenmp divide.c -o divide
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define BUFSIZE 65536 
#define FILENAME_LENGTH 1024 

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Format: ./divide [number_of_files] [output_file_prefix] [input_file_prefix] [number_to_divide_by]\n");
        exit(1);
    }

    const long number_of_files = atol(argv[1]);
    const unsigned int d = atoi(argv[4]);

    if (d == 0) {
        perror("The divisor cannot be 0\n");
        exit(1);
    }

    //size_t common_filesize = 0;

    // printf("Max threads: %d\n", omp_get_max_threads());

    #pragma omp parallel for
    for (int i = 0; i < number_of_files; i++) {
        // printf("Thread %d\n", omp_get_thread_num());
        char filename[FILENAME_LENGTH];
        size_t filesize;

        FILE *fp_input;
        FILE *fp_output;

        snprintf(filename, sizeof(filename), "%s%d", argv[3], i);
        fp_input = fopen(filename, "r");
        if (fp_input == NULL) {
            perror(filename);
            exit(1);
        }
        fseek(fp_input, 0, SEEK_END);
        filesize = ftell(fp_input);

        if (filesize == 0) {
            perror("File size must exceed 0\n");
            exit(1);
        }

        if (filesize % (BUFSIZE * sizeof(unsigned int)) != 0) {
            perror("Each file size must be divisible by BUFSIZE * sizeof(unsigned int)\n");
            exit(1);
        }

        /*if (common_filesize == 0) {
            common_filesize = filesize;
        } else {
            if (filesize != common_filesize) {
                perror("All files must have exactly the same size\n");
                exit(1);
            }
        }*/

        fseek(fp_input, 0, SEEK_SET);

        // Opening the file where we will write the result
        snprintf(filename, sizeof(filename), "%s%d", argv[2], i);
        fp_output = fopen(filename, "w");
        if (fp_output == NULL) {
            perror(filename);
            exit(1);
        }

        unsigned int buf[BUFSIZE];

        while (fread(buf, sizeof(unsigned int), BUFSIZE, fp_input) == BUFSIZE) {
            for (int b = 0; b < BUFSIZE; b++) {
                buf[b] /= d;
            }
            fwrite(buf, sizeof(unsigned int), BUFSIZE, fp_output);
        }

        fclose(fp_input);
        fclose(fp_output);
    }

    return 0;
}
