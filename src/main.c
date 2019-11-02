#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../inc/vptree.h"

#define POINTS 10000000
#define DIMENSIONS 10

void benchmark(){
    // Create a random X array
    double *input;
    double avg = 0;
    clock_t start, end;
    vptree *tree;


    printf("Points, Dimensions, Average time\n");
    for(int points = 100; points < 10000000; points = points * 10){

        for(int dimensions = 10; dimensions < 1000; dimensions = dimensions * 10){

            // We don't have enough RAM for that
            if(dimensions * points > 100000000)
                break;

            avg = 0;

            // Run tests
            for(int i=0; i < 10; i++){

                // Create input
                input = calloc(points * dimensions, sizeof(double));

                for(int j = 0; j < points * dimensions; j++){
                        input[j] = rand();
                }

                // Call and benchmark
                start = clock();
                buildvp(input, points, dimensions);

                end = clock();
                avg += ((double) (end - start)) / CLOCKS_PER_SEC;

            }

            // Free input - no need it's freed by the subroutine
            //free(input);

            avg = avg / 10;

            // Print results
            printf("%d, %d, %f\n", points, dimensions, avg);


        }

    }



}

int main()
{

    // Intialize random number generator
    time_t t;
    srand((unsigned) time(&t));

    benchmark();


    // Create a random X array
    /*
    double *X = calloc(POINTS * DIMENSIONS, sizeof(double));
    for(int i = 0; i < POINTS * DIMENSIONS; i++)
        X[i] = rand();

    clock_t start, end;
    start = clock();
    vptree *tree = buildvp(X, POINTS, DIMENSIONS);
    end = clock();
    printf("Root median: %f", tree->md);
    double elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nTime: %f", elapsed);*/

    //TODO: Add a function to visualize tree
    return 0;
}
