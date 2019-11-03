#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../inc/vptree.h"

// Function Prototypes
vptree * buildvp(double *X, int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double getMD(vptree * T);
double * getVP(vptree * T);
int getIDX(vptree * T);
vptree * build_tree(double *points, int *ids, int n, int d);
double * euclidean(double *point, double *points, int *ids, int n, int d);
void swap(double *a, double *b);
int partition (double arr[], int *ids, int low, int high);
double quickselect_median(double arr[], int *ids, int length);
double quickselect(double arr[], int *ids, int length, int idx);

// Application entry point
vptree * buildvp(double *X, int n, int d)
{
    // Allocate space for the index array
    int *ids = calloc(n, sizeof(int));

    // Build the initial ids array
    for (int i = 0; i < n; i++)
        ids[i] = i;

    // Call build_tree to get the pointer to the root of the tree
    vptree * node = build_tree(X, ids, n, d);

    // Free memory for ids and points arrays
    free(ids);
    return node;
}

// Function that recursively builds the binary tree and returns a pointer to its root
vptree * build_tree(double *points, int *ids, int n, int d)
{
    // Create node to be returned
    vptree *node = calloc(1, sizeof(vptree));

    // Check to end recursion: if points array is of size 0 - we are returning a leaf
    if (n == 1)
    {
        // Build node
        node->inner = NULL;
        node->outer = NULL;
        node->idx = ids[0];
        node->md = 0;
        node->vp = points + ids[0] * d;

        // Return node
        return node;
    }

    // Choose the last point in X as the vantage point
    int id = ids[n-1];
    double *point = (points + id*d);

    // Create array that holds euclidean distance of point from all other points
    double *distances = euclidean(point, points, ids, n-1, d);

    // At this point distances[i] indicates the distance of point with ids[i] from the vantage point
    // Find median and quicksort the distances (and the parallel ids array)
    double median = quickselect_median(distances, ids, n-1);

    // At this point distances are sorted. The ones left of the median are smaller than it and the ones on the right larger
    // Also the point with index ids[i] has distance distances[i] from the vantage point

    // Sort points into two new arrays
    // Calculate array sizes for subtrees. Values up to and equal to the median go on the inner tree
    int innerLength = 0;
    for (int i = 0; i < n-1; i++)
    {
        if(distances[i] <= median)
        {
            innerLength++;
        }
    }
    int outerLength = n - 1 - innerLength;

    // Get pointers to the parts of ids which correspond to points that belong on the inner and outer subtrees
    int *innerIDs = ids;
    int *outerIDs = ids + innerLength;

    // Set node fields
    node->md = median;
    node->vp = point;
    node->idx = id;

    // De-allocate unused memory
    free(distances);

    // Calculate subtrees
    if(innerLength > 0){
        node->inner = build_tree(points, innerIDs, innerLength, d);
    }
    else{
        node->inner = NULL;
    }

    if(outerLength > 0){
        node->outer = build_tree(points, outerIDs, outerLength, d);
    }
    else{
        node->outer = NULL;
    }

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
double * euclidean(double *point, double *points, int *ids, int n, int d){

    double *distances = (double*)calloc(n, sizeof(double));

    double accumulator = 0;

    for (int i = 0; i < n; i++){
        accumulator = 0;
        for (int j = 0; j < d; j++){
            accumulator += (point[j] - *(points + ids[i] * d + j)) * (point[j] - *(points + ids[i] * d + j));
        }
        distances[i] = sqrt(accumulator);
    }
    return distances;
}

// A utility function to swap two elements
void swap_double(double *a, double *b){
    double t = *a;
    *a = *b;
    *b = t;
}

// A utility function to swap two elements
void swap_int(int *a, int *b){
    int t = *a;
    *a = *b;
    *b = t;
}


// QuickSort Partition function. low and high define the range of indexes in arr where partition should work
int partition (double arr[], int *ids, int low, int high){

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
            // Swap distances and corresponding point ids
            swap_double(&arr[i], &arr[j]);
            swap_int(&ids[i], &ids[j]);
        }
    }

    // Finally place pivot in its correct position in the array and return the position as the middle point
    swap_double(&arr[i + 1], &arr[high]);
    swap_int(&ids[i + 1], &ids[high]);
    return (i + 1);
}

// Returns the median using the QuickSelect algorithm
double quickselect_median(double arr[], int *ids, int length){

    if (length % 2 == 1){
        return quickselect(arr, ids, length, (length+1)/2);
    }
    else{
        return 0.5 * (quickselect(arr, ids, length, length/2) + quickselect(arr, ids, length, length/2 + 1));
    }
}

// Returns the idx-th element of arr when arr is sorted
// idx is the index (starting from 1) of the point we want to find when the array is sorted. For the median idx should be the middle one (i.e (length+1)/2 for odd lengths etc)
double quickselect(double arr[], int *ids, int length, int idx){

    // Check to end recursion
    if (length == 1){
        return arr[0];
    }

    // Select last array element as pivot
    double pivot = arr[length - 1];
    // Get index of pivot after we partition the array
    int pivotIndex = partition(arr, ids, 0, length - 1);

    // Create the higher and lower arrays that occur after partitioning in QuickSort fashion
    int lowerLength = pivotIndex;
    pivotIndex++;
    int higherLength = (length - (lowerLength + 1));

    // At this point pivotIndex, lowerLength and higherLength all start from 1 not 0

    // Variable to store result of following recursive calls
    double result = 0;

    // This means that the point we're looking (median in our case) is in the lower partition
    if (idx <= lowerLength){
        result = quickselect(arr, ids, lowerLength, idx);
    }
    // This means that the median is our pivot point
    else if(idx == pivotIndex){
        result = pivot;
    }
    // This means that the median is in the higher partition
    else{
        result = quickselect(arr + pivotIndex, ids + pivotIndex, higherLength, idx - pivotIndex);
    }

    // Return result
    return result;
}
