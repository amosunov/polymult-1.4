/* This program takes as an input prefixes of input binary files, sums their contents and saves them into an output file

EXAMPLE: ./sum 32 ~/res ~/a ~/b

First input files have names a0, a1, ..., a31
Second input files have names b0, b1, ..., b31
The resulting files are called res0, res1, ..., res31
If we view each file as a huge vector, then res0=a0+b0, res1=a1+b1, ..., res31=a31+b31 

Notice that all files must store binary unsigned int values, have exactly the same nonzero size, and this size must be divisible by BUFSIZE

Compile with clang -fopenmp sum.c -o sum
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#define BUFSIZE 65536 
#define FILENAME_LENGTH 1024 

int main(int argc, char *argv[]) {
    #ifdef WITH_PARI
    pari_init(8000000, 500000);
    #endif

    if (argc < 5) {
        printf("Format: ./sum [number_of_files] [output_file_prefix] [input_file_prefix1] [input_file_prefix2] ...\n");
        exit(1);
    }

    const long number_of_files = atol(argv[1]);

    //size_t common_filesize = 0;

    // printf("Max threads: %d\n", omp_get_max_threads());

    #pragma omp parallel for
    for (int i = 0; i < number_of_files; i++) {
        // printf("Thread %d\n", omp_get_thread_num());
        char filename[FILENAME_LENGTH];
        size_t filesize;

        FILE *fp[argc - 2];

        int j;
        for (j = 1; j < argc - 2; j++) {
            snprintf(filename, sizeof(filename), "%s%d", argv[j + 2], i);
            fp[j] = fopen(filename, "r");
            if (fp[j] == NULL) {
                perror(filename);
                exit(1);
            }
            fseek(fp[j], 0, SEEK_END);
            filesize = ftell(fp[j]);

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

            fseek(fp[j], 0, SEEK_SET);
        }

        // Opening the file where we will write the result
        snprintf(filename, sizeof(filename), "%s%d", argv[2], i);
        fp[0] = fopen(filename, "w");
        if (fp[0] == NULL) {
            perror(filename);
            exit(1);
        }

        unsigned int buf[BUFSIZE];
        unsigned int result[BUFSIZE];

        while (fread(result, sizeof(unsigned int), BUFSIZE, fp[1]) == BUFSIZE) {
           for (j = 2; j < argc - 2; j++) {
                fread(buf, sizeof(unsigned int), BUFSIZE, fp[j]);
                for (int b = 0; b < BUFSIZE; b++) {
                    result[b] += buf[b];
                }
           }
           fwrite(result, sizeof(unsigned int), BUFSIZE, fp[0]);
        }

        /*
        int k = 0;

        while (fread(&value1, sizeof(unsigned int), 1, fp1) == 1 && fread(&value2, sizeof(unsigned int), 1, fp2)) {
            //printf("h(-%u)=%u\n", 120*k + 23, value1+value2);
            #ifdef WITH_PARI
            pari_sm av = avma;
            GEN D = stoi(-(120*k + 23));
            if (Z_issquarefree(D)) {
                GEN class_number = gel(quadclassunit0(D, 0, NULL, DEFAULTPREC), 1);
                if ((value1 + value2) != itou(class_number)) {
                 printf("k=%d, D=%d, CALCULATED: %u, ACTUAL: %u\n", k, -(120*k + 23), value1 + value2, (unsigned int) itou(class_number));
                }
            }
            avma = av;
            #endif
            k++;
        }*/

        fclose(fp[0]);
        for (j = 1; j < argc - 2; j++) {
            fclose(fp[j]);

            #ifdef DELETE_INPUT_FILES
            snprintf(filename, sizeof(filename), "%s%d", argv[j + 2], i);
            remove(filename);
            #endif
        }
    }

    #ifdef WITH_PARI
    pari_close();
    #endif

    return 0;
}
