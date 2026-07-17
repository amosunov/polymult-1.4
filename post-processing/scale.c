/*
Takes as input a prefix of input binary files, divides or multiplies their contents (unsigned integers) by a fixed positive unsigned integer, and saves the resulting values into an output file

EXAMPLE: ./sum 32 ~/res ~/a divide 3

The input files have names a0, a1, ..., a31
The output files are called res0, res1, ..., res31
If we view both files as huge vectors, then res0=(1/3)*a0, res1(1/3)*a1, ..., res31=(1/3)*a31 

Notice that all files must store binary unsigned int values, have exactly the same nonzero size, and this size must be divisible by BUFSIZE

Compile with clang -fopenmp scale.c -o scale
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define BUFSIZE 65536 
#define FILENAME_LENGTH 1024 

#define MULTIPLY 0 
#define DIVIDE 1

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Format: ./divide [number_of_files] [output_file_prefix] [input_file_prefix] [operation] [number_to_multiply_or_divide_by]\n");
        exit(1);
    }

    const long number_of_files = atol(argv[1]);
    const int operation = strcmp(argv[4], "multiply");
    const unsigned int d = atoi(argv[5]);

    if (d == 0) {
        perror("The divisor cannot be 0\n");
        exit(1);
    }

    const int same_input_output_files = (strcmp(argv[2], argv[3]) == 0);

    #pragma omp parallel for
    for (int i = 0; i < number_of_files; i++) {
        char filename[FILENAME_LENGTH];
        size_t filesize;
        unsigned int buf[BUFSIZE];
        size_t total_read;
        int b;

        if (same_input_output_files) {
            FILE *fp;
            snprintf(filename, sizeof(filename), "%s%d", argv[3], i);
            fp = fopen(filename, "r+b");
            if (fp == NULL) {
                perror(filename);
                exit(1);
            }

            total_read = fread(buf, sizeof(unsigned int), BUFSIZE, fp);
            while (total_read) {
                if (operation == DIVIDE) {
                    for (b = 0; b < total_read; b++) {
                        buf[b] /= d;
                    }
                } else { // otherwise multiply
                    for (b = 0; b < total_read; b++) {
                        buf[b] *= d;
                    }
                }
                fseek(fp, -(long) (total_read * sizeof(unsigned int)), SEEK_CUR);
                fwrite(buf, sizeof(unsigned int), total_read, fp);
                total_read = fread(buf, sizeof(unsigned int), BUFSIZE, fp);
            }

            fclose(fp);
        } else {
            FILE *fp_input;
            FILE *fp_output;

            snprintf(filename, sizeof(filename), "%s%d", argv[3], i);
            fp_input = fopen(filename, "r");
            if (fp_input == NULL) {
                perror(filename);
                exit(1);
            }

            // Opening the file where we will write the result
            snprintf(filename, sizeof(filename), "%s%d", argv[2], i);
            fp_output = fopen(filename, "w");
            if (fp_output == NULL) {
                perror(filename);
                exit(1);
            }

            total_read = fread(buf, sizeof(unsigned int), BUFSIZE, fp_input);
            while (total_read) {
                if (operation == DIVIDE) {
                    for (b = 0; b < total_read; b++) {
                        buf[b] /= d;
                    }
                } else { // otherwise multiply
                    for (b = 0; b < total_read; b++) {
                        buf[b] *= d;
                    }
                }
                fwrite(buf, sizeof(unsigned int), total_read, fp_output);
                total_read = fread(buf, sizeof(unsigned int), BUFSIZE, fp_input);
            }

            fclose(fp_input);
            fclose(fp_output);
        
            #ifdef DELETE_INPUT_FILES
            snprintf(filename, sizeof(filename), "%s%d", argv[3], i);
            remove(filename);
            #endif
        }
    }

    return 0;
}
