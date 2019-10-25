#include <stdio.h>
#include <stdlib.h>
#include "../inc/vptree.h"

#define POINTS 100000
#define DIMENSIONS 1000

int main()
{

    // Intialize random number generator
    time_t t;
    srand((unsigned) time(&t));

    // Create a random X array
    double *X = calloc(POINTS * DIMENSIONS, sizeof(double));
    for(int i = 0; i < POINTS * DIMENSIONS; i++)
        X[i] = rand();
    vptree *tree = buildvp(X, POINTS, DIMENSIONS);
    printf("Root median: %f", tree->md);

    //TODO: Add a function to visualize tree
    return 0;
}
