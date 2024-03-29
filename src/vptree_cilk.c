#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "../inc/vptree.h"

// Number of threads for Cilk to use
#define THREADS "8"
// Threshold of points to switch to sequential execution. I recommend turning this off (setting to 0) and let cilk manage it
#define POINT_THRESHOLD 0
// Development flags to switch execution mode from serial to parallel for distance calculation and subtree creation
#define PARALLELDIS true
#define PARALLELSUB true

// Function Prototypes
vptree * buildvp(double *X, int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double getMD(vptree * T);
double * getVP(vptree * T);
int getIDX(vptree * T);
vptree * build_tree(double *points, int *ids, int n, int d);
void euclidean(double *point, double *points, double *distances, int n, int d);
void swap(double *a, double *b);
int partition (double arr[], int low, int high);
double quickselect_median(double arr[], int length);
double quickselect(double arr[], int length, int idx);

// Application entry point
vptree * buildvp(double *X, int n, int d)
{
    // Allocate space for the index array
    int *ids = calloc(n, sizeof(int));

    // Build the initial ids array
    for (int i = 0; i < n; i++)
        ids[i] = i;

    // Call build_tree to get the pointer to the root of the tree
    return build_tree(X, ids, n, d);
}

// Function that recursively builds the binary tree and returns a pointer to its root
vptree * build_tree(double *points, int *ids, int n, int d)
{
	// Set number of threads. I recommend commenting this out and leaving Cilk to decide this
	//__cilkrts_set_param("nworkers", THREADS_MAX);

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
        node->vp = calloc(d, sizeof(double));
        memcpy(node->vp, points, sizeof(double) * d);

        // Free memory for ids and points arrays
        free(ids);
        free(points);

        // Return node
        return node;
    }

    // Choose the last point in X as the vantage point
    double *point = (points + (n-1)*d);
    double id = ids[n-1];

    // Create array that holds euclidean distance of point from all other points
    double *distances = calloc(n-1, sizeof(double));

	// Calculate distances in parallel if possible using work-construct, else do it sequentially (logic in euclidean())
	euclidean(point, points, distances, n-1, d);

    // At this point distances[i] indicates the distance of point i in points from the vantage point
    // Find median by creating a copy of distances and passing it to QuickSelect
    double *distancesCopy = calloc(n-1, sizeof(double));
    memcpy(distancesCopy, distances, sizeof(double) * (n-1));
    double median = quickselect_median(distancesCopy, n-1);
    free(distancesCopy);

    // Sort points into two new arrays
    // Calculate array sizes for subtrees. Values up to and equal to the median go on the inner tree
    int innerLength = 0;
    for (int i = 0; i < n-1; i++){
        if(distances[i] <= median){
            innerLength++;
        }
    }
    int outerLength = n - 1 - innerLength;

    // Pointers to keep track of inner and outer arrays content while sorting points
    int innerPointer = 0;
    int outerPointer = 0;

    // Create arrays for inner and outer points. Arrays contain the points and a list of ids (one for each point)
    double *innerPoints = calloc(innerLength * d, sizeof(double));
    double *outerPoints = calloc(outerLength * d, sizeof(double));
    int *innerIDs = calloc(innerLength, sizeof(int));
    int *outerIDs = calloc(outerLength, sizeof(int));

    // Sort points to inner and outer subtree
    for (int i = 0; i < n-1; i++){
        if(distances[i] <= median){
            memcpy(innerPoints + innerPointer * d, points + i*d, sizeof(double) * d);
            innerIDs[innerPointer] = ids[i];
            innerPointer++;
        }
        else{
            memcpy(outerPoints + outerPointer * d, points + i*d, sizeof(double) * d);
            outerIDs[outerPointer] = ids[i];
            outerPointer++;
        }
    }

    // Set node fields
    // Copy the point into vp because we will call free(points) that will also free(point)
    node->vp = calloc(d, sizeof(double));
    node->md = median;
    memcpy(node->vp, point, sizeof(double) * d);
    node->idx = id;

    // De-allocate unused memory
    free(points);
    free(distances);
    free(ids);

    // Calculate subtrees in parallel if possible
    if((innerLength > 0) && (outerLength > 0)){
        node->outer = cilk_spawn build_tree(outerPoints, outerIDs, outerLength, d);
        node->inner = build_tree(innerPoints, innerIDs, innerLength, d);
        cilk_sync;
    }
    else if ((innerLength > 0) && !(outerLength > 0)){
    	node->inner = build_tree(innerPoints, innerIDs,innerLength, d);
    	node->outer = NULL;
    }
    else if (!(innerLength > 0) && (outerLength > 0)){
    	node->outer = build_tree(outerPoints, outerIDs, outerLength, d);
		node->inner = NULL;
	}
    else{
        node->outer = NULL;
        node->inner = NULL;
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

// Calculates the distances of all points from point and writes them to an array. If possible use Cilk Parallel Array operations
void euclidean(double *point, double *points, double *distances, int n, int d)
{
	// Accumulator array for parallel execution
    double *accumulator;
    // Scalar accumulator for sequential execution
    int accumulator_seq = 0;

    // Decide if point calculation should happen in parallel or not
    if((n-1 > POINT_THRESHOLD) && (PARALLELDIS == true))
    {
    	accumulator = calloc(n * d, sizeof(double));

		cilk_for (int i = 0; i < n; i++)
		{
			// Vectorized calculations
			accumulator[i*d:d] = point[0:d] - points[i*d:d];
			accumulator[i*d:d] = accumulator[i*d:d] * accumulator[i*d:d];
			// Use a reducer to sum up all elements of the accumulator vector and get its sqrt
			distances[i] = sqrt(__sec_reduce_add (accumulator[i*d:d]));
		}
	    free(accumulator);

    }else{
        for (int i = 0; i < n; i++)
        {
            accumulator_seq = 0;
            for (int j = 0; j < d; j++)
            {
                accumulator_seq += (point[j] - *(points + i * d + j)) * (point[j] - *(points + i * d + j));
            }
            //TODO: If the problem gets too big in size and accumulator values rise we will have to perform the sqrt operation
            // on every iteration to avoid overflowing
            distances[i] = sqrt(accumulator_seq);
        }
    }
    return;
}

// A utility function to swap two elements
void swap(double *a, double *b){
    double t = *a;
    *a = *b;
    *b = t;
}

// QuickSort Partition function. low and high are the range of indexes in arr where partition should work
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
    pivotIndex++;
    int higherLength = (length - (lowerLength + 1));

    // At this point pivotIndex, lowerLength and higherLength all start from 1 not 0
    double *lower = calloc(lowerLength, sizeof(double));
    double *higher = calloc(higherLength, sizeof(double));
    memcpy(lower, arr, sizeof(double) * lowerLength);
    memcpy(higher, arr + pivotIndex, sizeof(double) * higherLength);

    // Variable to store result of following recursive calls
    double result = 0;

    // This means that the point we're looking (median in our case) is in the lower partition
    if (idx <= lowerLength){
        result = quickselect(lower, lowerLength, idx);
    }
    // This means that the median is our pivot point
    else if(idx == pivotIndex){
        result = pivot;
    }
    // This means that the median is in the higher partition
    else{
        result = quickselect(higher, higherLength, idx - pivotIndex);
    }

    // Free memory allocated to lower and higher
    free(lower);
    free(higher);

    // Return result
    return result;
}
