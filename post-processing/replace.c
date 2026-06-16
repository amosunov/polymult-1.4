// This program takes as an input an unsigned integer k, and index n, and a binary file containing unsigned integers, and replaces the n-th entry in the binary file with k.
// Compile with clang replace.c -o replace

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: ./replace [filename] [index] [new_value]\n");
        return 1;
    }

    const char *filename = argv[1];
    long n = atol(argv[2]);
    unsigned int k = atoi(argv[3]);

    FILE *fp = fopen(filename, "r+b");
    if (fp == NULL) {
        perror("fopen");
        return 1;
    }

    /* Move to the n-th unsigned integer */
    if (fseek(fp, n * sizeof(unsigned int), SEEK_SET) != 0) {
        perror("fseek");
        fclose(fp);
        return 1;
    }

    /* Overwrite the value */
    if (fwrite(&k, sizeof(unsigned int), 1, fp) != 1) {
        perror("Unable to write\n");
        fclose(fp);
        return 1;
    }

    fclose(fp);
    return 0;
}
