/*
Merges data containing class numbers for fundamental discriminants with |D| = a1, a2, ..., ak (mod m_{big}) into a single file with |D| = a_1 (mod m_{small}). In order for the program to work properly, must have 0 <= a1 < m_{small}$, from which a1, a2=a1+m_{small}, a3=a1+2*m_{small} are derived. If a file for a certain congruence class is missing, zero is recorded in place of the actual class number.

In order for the program to work properly all existing files must have exactly the same size, and this size must be divisible by BUFSIZE.

Compile with clang -fopenmp merge.c -o merge
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define BUFSIZE 65536
#define FILENAME_LENGTH 1024

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Format: ./merge [path] [number_of_files] [a] [small_m] [big_m]\n");
        exit(1);
    }

    char *path = argv[1];
    const int number_of_files = atoi(argv[2]);
    const int a = atoi(argv[3]);
    const int small_m = atoi(argv[4]);
    const int big_m = atoi(argv[5]);

    if (big_m % small_m != 0) {
        perror("small_m must divide big_m\n");
        exit(1);    
    }

    if (a < 0 || a >= small_m) {
        perror("a must satisfy 0 <= a < small_m\n");
        exit(1);
    }

    const int total_congruence_classes = big_m/small_m;

    long number_of_discs_per_file = 0;

    #pragma omp parallel for
    for (int i = 0; i < number_of_files; i++) {
        char filename[FILENAME_LENGTH];
        FILE *fp[total_congruence_classes + 1];

        int a_var = a, j;
        long number_of_discs_per_file = 0;
        for (j = 1; j <= total_congruence_classes; j++) {
            snprintf(filename, sizeof(filename), "%s/h%dmod%d.%d", path, a_var, big_m, i);
            #ifdef DEBUG
            printf("Opening the file %s\n", filename);
            #endif
            fp[j] = fopen(filename, "r");
            a_var += small_m;
            if (fp[j] != NULL) {
                #pragma omp critical
                {
                    if (number_of_discs_per_file == 0) {
                        fseek(fp[j], 0, SEEK_END);
                        number_of_discs_per_file = ftell(fp[j]);
                        fseek(fp[j], 0, SEEK_SET);
                    }
                }
            }
        }
        number_of_discs_per_file /= sizeof(unsigned int);

        #ifdef DEBUG
        printf("There are %ld discriminants per file\n", number_of_discs_per_file);
        #endif

        // Opening the file where we will write the result
        snprintf(filename, sizeof(filename), "%s/h%dmod%d.%d", path, a, small_m, i);
        fp[0] = fopen(filename, "w");

        unsigned int buf[total_congruence_classes + 1][BUFSIZE];
        unsigned int result[total_congruence_classes * BUFSIZE];

        size_t total_read;

        #ifdef DEBUG
        unsigned long D = small_m * total_congruence_classes * number_of_discs_per_file * i + a;
        printf("For the file #%d the first discriminant is %ld\n", i, D);
        #endif

        for (j = 1; j <= total_congruence_classes; j++) {
            if (fp[j] != NULL) {
                memset(buf[j], 0, BUFSIZE * sizeof(unsigned int));
            }
        }

        while (1) {
            total_read = 0;
            for (j = 1; j <= total_congruence_classes; j++) {
                if (fp[j] != NULL) {
                    total_read += fread(buf[j], sizeof(unsigned int), BUFSIZE, fp[j]);
                } else {
                    memset(buf[j], 0, BUFSIZE * sizeof(unsigned int));
                }
            }

            if (total_read == 0) break;

            memset(result, 0, total_congruence_classes * BUFSIZE * sizeof(unsigned int));

            int r = 0;
            for (int k = 0; k < BUFSIZE; k++) {
                for (j = 1; j <= total_congruence_classes; j++) {
                    #ifdef DEBUG
                    printf("H(%d)=%u\n", D, buf[j][k]);
                    D += small_m;
                    #endif
                    result[r] = buf[j][k];
                    r++;
                }
            }

            fwrite(result, sizeof(unsigned int), total_congruence_classes * BUFSIZE, fp[0]);
        }

        for (j = 0; j <= total_congruence_classes; j++) {
            fclose(fp[j]);
        }
    }
}
