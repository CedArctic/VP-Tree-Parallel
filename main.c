#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../inc/vptree.h"

#define POINTS 10000
#define DIMENSIONS 700

// Function Prototypes
vptree * buildvp(double *X, int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double getMD(vptree * T);
double * getVP(vptree * T);
int getIDX(vptree * T);
double * euclidean(double *point, double *points, int n, int d);
void swap(double *a, double *b);
int partition (double arr[], int low, int high);
double quickselect_median(double arr[], int length);
double quickselect(double arr[], int length, int idx);


// Build the tree
vptree * buildvp(double *X, int n, int d){

    // Create node to be returned
    vptree *node = malloc(sizeof(vptree));

    // Check to end recursion: if points array is of size 0 - we are returning a leaf
    if (n == 1){
        node->inner = NULL;
        node->outer = NULL;
        node->md = 0;
        node->vp = calloc(d, sizeof(double));
        memcpy(node->vp, X, sizeof(double) * d);
        return node;
    }

    // Choose a random point as the vantage point: in this case we take the last one in X
    double *point = calloc(d, sizeof(double));
    // Copy the point from the original matrix to a new vector
    memcpy(point, (X + (n-1)*d), sizeof(double) * d);

    // Copy all other points to a new vector (i.e all of X from rows 1 to n-1)
    double *points = calloc((n-1) * d, sizeof(double));
    memcpy(points, X, sizeof(double) * d * (n-1));

    // Create array that holds euclidean distance of point from all other points
    double *distances = euclidean(point, points, n-1, d);


    // At this point distances[i] indicates the distance of point i in points from the vantage point
    // Find median by creating a copy of distances and passing it to QuickSelect
    double *distancesCopy = calloc(n-1, sizeof(double));
    memcpy(distancesCopy, distances, sizeof(double) * (n-1));
    double median = quickselect_median(distancesCopy, n-1);

    // Sort points into two new arrays
    // Calculate array sizes for subtrees. In the event that n-1 is an odd number, the point with median distance goes to the outer subtree
    int innerLength = (int)floor((n-1) / 2);
    int outerLength = n - 1 - innerLength;
    // ERROR: Outer length calculation seems to be off when using ceil(). That's why above workaround is used
    //int outerLength = (int)ceil((n-1) / 2);

    // Pointers to keep track of inner and outer arrays content while sorting points
    int innerPointer = 0;
    int outerPointer = 0;

    // Create and return node
    double *innerPoints = calloc(innerLength * d, sizeof(double));
    double *outerPoints = calloc(outerLength * d, sizeof(double));

    // Sort points
    for (int i = 0; i < n-1; i++){
        if(distances[i] < median){
            memcpy(innerPoints + innerPointer * d, X + i*d, sizeof(double) * d);
            innerPointer++;
        }
        else{
            memcpy(outerPoints + outerPointer * d, X + i*d, sizeof(double) * d);
            outerPointer++;
        }
    }

    //TODO: Maybe assert that innerPointer == innerLength - 1 at this point

    if(innerLength > 0){
       node->inner = buildvp(innerPoints, innerLength, d);
    }
    else{
        node->inner = NULL;
    }

    if(outerLength > 0){
       node->outer = buildvp(outerPoints, outerLength, d);
    }
    else{
        node->outer = NULL;
    }
    node->md = median;
    node->vp = point;

    // De-allocate unused memory
    free(points);
    free(distances);
    free(distancesCopy);
    free(innerPoints);
    free(outerPoints);

    return node;
}

// Return vantage-point subtree with points inside radius
vptree * getInner(vptree * T){
    return T->inner;
}

// Return vantage-point subtree with points outside radius
vptree * getOuter(vptree * T){
    return T->outer;
}

// Return median of distances to vantage point
double getMD(vptree * T){
    return T->md;
}

// Return the coordinates of the vantage point
double * getVP(vptree * T){
    return T->vp;
}

// Return the index of the vantage point
int getIDX(vptree * T){
    return T->idx;
}

// Returns pointer to an array that contains the distances of all points from point
double * euclidean(double *point, double *points, int n, int d){

    double *distances = (double*)calloc(n, sizeof(double));

    double accumulator = 0;

    for (int i = 0; i < n; i++){
        accumulator = 0;
        for (int j = 0; j < d; j++){
            //printf("Current (%d) points element: %.2f\n", i*d+j, *(points + i * d + j));
            accumulator += pow((point[j] - *(points + i * d + j)), 2);
        }
        distances[i] = sqrt(accumulator);
    }
    return distances;
}

// A utility function to swap two elements
void swap(double *a, double *b){
    double t = *a;
    *a = *b;
    *b = t;
}

// QuickSort Partition function
// low and high are the range of indexes in arr where partition should work
int partition (double arr[], int low, int high){

    // Select a pivot and initialize flag to position of smallest element before pivot
    double pivot = arr[high];
    int i = (low - 1);

    // Go through the array examining each element
    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot, increment i and swap it out with the one currently pointed by i
        if (arr[j] < pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }

    // Finally place pivot in its correct position in the array and return the position as the middle point
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// Returns the median using the QuickSelect algorithm
double quickselect_median(double arr[], int length){

    if (length % 2 == 1){
        return quickselect(arr, length, (length+1)/2);
    }
    else{
        return 0.5 * (quickselect(arr, length, length/2) + quickselect(arr, length, length/2 + 1));
    }
}

// Returns the idx-th element of arr when arr is sorted
// idx is the index (starting from 1) of the point we want to find when the array is sorted. For the median idx should be the middle one (i.e (length+1)/2 for odd lengths etc)
double quickselect(double arr[], int length, int idx){

    // Check to end recursion
    if (length == 1){
        return arr[0];
    }

    // Select last array element as pivot
    double pivot = arr[length - 1];
    // Get index of pivot after we partition the array
    int pivotIndex = partition(arr, 0, length - 1);

    // Create the higher and lower arrays that occur after partitioning in QuickSort fashion
    int lowerLength = pivotIndex;
    int higherLength = (length - (pivotIndex + 1));
    double *lower = calloc(lowerLength, sizeof(double));
    double *higher = calloc(higherLength, sizeof(double));
    memcpy(lower, arr, sizeof(double) * lowerLength);
    memcpy(higher, arr + pivotIndex + 1, sizeof(double) * higherLength);

    // Variable to store result of following recursive calls
    double result = 0;

    // This means that the point we're looking (median in our case) is in the lower partition
    if (idx <= lowerLength){
        result = quickselect(lower, lowerLength, idx);
    }
    // This means that the median is our pivot point
    else if(idx <= lowerLength + 1){
        result = pivot;
    }
    // This means that the median is in the higher partition
    else{
        result =  quickselect(higher, higherLength, idx - lowerLength - 1);
    }

    // Free memory allocated to lower and higher
    free(lower);
    free(higher);

    // Return result
    return result;
}

int main()
{

    // Intialize random number generator
    time_t t;
    srand((unsigned) time(&t));

    // Create a random X array
    double *X = calloc(POINTS * DIMENSIONS, sizeof(double));
    for(int i = 0; i < POINTS * DIMENSIONS; i++)
        X[i] = rand() % 50;
    vptree *tree = buildvp(X, POINTS, DIMENSIONS);
    printf("Root median: %f", tree->md);
    


    /*
    double arr[10] = {12,2,3,3,5,19,7,8,9,10};
    vptree *tree = buildvp(arr, 5, 2);
    printf("Root median: %f", tree->md);
    */

    //TODO: Add a function to visualize tree
    return 0;
}
