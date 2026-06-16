/*
Uses PARI/GP to check if class numbers h(D) for |D|=a (mod m) in a given file were computed correctly
clang check.c -I/opt/homebrew/include -L/opt/homebrew/lib -lpari -o check

In order for the program to work correcrly, file size must be divisible by BUFSIZE
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <pari/pari.h>

#define BUFSIZE 65536 

int main(int argc, char *argv[]) {
    pari_init(8000000, 500000);

    if (argc != 5) {
        printf("Format: ./check [file_with_class_numbers] [first_discriminant] [a] [m]\n");
        exit(1);
    }

    const char * file_with_class_numbers = argv[1];
    long D = atol(argv[2]);
    const long a = atol(argv[3]);
    const long m = atol(argv[4]);

    if ((D - a) % m != 0) {
        perror("First discriminant must be congruent to a modulo m\n");
        exit(1);
    }

    // To check if the discriminant is fundanemtal we need to divide D=0 (mod 4) by 4
    long r = (D % 4 == 0) ? 4 : 1;

    FILE *fp;

    fp = fopen(file_with_class_numbers, "r");
    if (fp == NULL) {
        perror(file_with_class_numbers);
        exit(1);
    }

    unsigned int h[BUFSIZE];

    while (fread(h, sizeof(unsigned int), BUFSIZE, fp) == BUFSIZE) {
        for (int j = 0; j < BUFSIZE; j++) {
            pari_sp av = avma;
            if ((D % 120) != 71 && (D % 120) != 119) {
                GEN D_pari = stoi((-1)*D);
                GEN D_divide_by_r = stoi(D/r);
                if (Z_issquarefree(D_divide_by_r)) {
                    GEN class_number = gel(quadclassunit0(D_pari, 0, NULL, DEFAULTPREC), 1);
                    if (h[j] != itou(class_number)) {
                        printf("\tERROR: h(%ld)=%ld, but the file has %d\n", (-1)*D, itou(class_number), h[j]);
                    }
                }
            } else {
                if (h[j] != 0) {
                    printf("\tERROR: The file has %u, but for D=%lu (mod 120) must have 0\n", h[j], D % 120);
                }
            }
            D += m;
            avma = av;
        }
    }

    fclose(fp);

    pari_close();

    return 0;
}
