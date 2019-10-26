#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "../inc/vptree.h"

// Number of threads in parallel distance calculation
#define THREADS 8
// Threshold of points to switch to sequential execution
#define POINT_THRESHOLD 100
// Threshold of maximum live threads simultaneously
#define THREADS_MAX 100
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
void euclidean(double *point, double *points, double *distances, int n, int d);
void swap(double *a, double *b);
int partition (double arr[], int low, int high);
double quickselect_median(double arr[], int length);
double quickselect(double arr[], int length, int idx);
void *buildvp_wrapper(void *arg);

// Flag used to detect if build_vp has already been called. If it has not, X is the original input array.
// If it has it means that X is the points vector with an idx vector extended to it at the end
bool runFlag = false;

// Counter to keep track of live threads
int threadCount = 0;

// Thread-safe function to alter the live thread count
void modThreadCount(int n){
    #pragma omp atomic
    threadCount += n;
}

// Function that recursively builds the binary tree
vptree * buildvp(double *X, int n, int d)
{
    // Enable OpenMP Nested Parallelism
    omp_set_nested(true);

    // Allocate space for the index array
    double *ids = calloc(n, sizeof(double));

    // If runFlag is true -> not first execution -> X has an ids array at its end
    // Else if we're on the first run, we're going to generate the ids array
    if (runFlag == true)
    {
        memcpy(ids, X + n * d, sizeof(double) * n);
    }
    else
    {
        for (int i = 0; i < n; i++)
            ids[i] = i;
    }

    // Set run flag to true
    runFlag = true;

    // Create node to be returned
    vptree *node = calloc(1, sizeof(vptree));

    // Check to end recursion: if points array is of size 0 - we are returning a leaf
    if (n == 1)
    {
        node->inner = NULL;
        node->outer = NULL;
        node->idx = ids[0];
        node->md = 0;
        node->vp = calloc(d, sizeof(double));
        memcpy(node->vp, X, sizeof(double) * d);
        return node;
    }

    // Choose the last point in X as the vantage point
    double *point = calloc(d, sizeof(double));
    double id = ids[n-1];

    // Copy the point from the original matrix to a new vector
    memcpy(point, (X + (n-1)*d), sizeof(double) * d);

    // Copy all other points to a new vector (i.e all of X from rows 1 to n-1)
    double *points = calloc((n-1) * d, sizeof(double));
    memcpy(points, X, sizeof(double) * d * (n-1));

    // Create array that holds euclidean distance of point from all other points
    double *distances = calloc(n-1, sizeof(double));

    // Calculate distances in parallel if possible using work-construct, else do it sequentially (logic in euclidean())
    euclidean(point, points, distances, n-1, d);

    // At this point distances[i] indicates the distance of point i in points from the vantage point
    // Find median by creating a copy of distances and passing it to QuickSelect
    double *distancesCopy = calloc(n-1, sizeof(double));
    memcpy(distancesCopy, distances, sizeof(double) * (n-1));
    double median = quickselect_median(distancesCopy, n-1);

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
    //TODO: Perhaps use distancesCopy to reduce the above linear scan to half

    // Pointers to keep track of inner and outer arrays content while sorting points
    int innerPointer = 0;
    int outerPointer = 0;

    // Create arrays for inner and outer points. Arrays contain the points and a list of ids (one for each point)
    double *innerPoints = calloc(innerLength * d + innerLength, sizeof(double));
    double *outerPoints = calloc(outerLength * d + outerLength, sizeof(double));

    // Sort points
    for (int i = 0; i < n-1; i++)
    {
        if(distances[i] <= median)
        {
            memcpy(innerPoints + innerPointer * d, X + i*d, sizeof(double) * d);
            innerPoints[innerLength * d + innerPointer] = ids[i];
            innerPointer++;
        }
        else
        {
            memcpy(outerPoints + outerPointer * d, X + i*d, sizeof(double) * d);
            outerPoints[outerLength * d + outerPointer] = ids[i];
            outerPointer++;
        }
    }

    // Build subtrees in parallel or sequentially
    if((PARALLELSUB == true) && (THREADS_MAX - threadCount >= 2))
    {
        modThreadCount(2);
        #pragma omp parallel shared(node, d, innerPoints, innerLength, outerPoints, outerLength)
        {
            #pragma omp sections nowait
            {
                // Create threads
                #pragma omp section
                if(innerLength > 0)
                {
                    node->inner = buildvp(innerPoints, innerLength, d);
                }

                #pragma omp section
                if(outerLength > 0)
                {
                    node->outer = buildvp(outerPoints, outerLength, d);
                }
            }
        }
        modThreadCount(-2);
    }
    else
    {
        if(innerLength > 0)
        {
            node->inner = buildvp(innerPoints, innerLength, d);
        }
        if(outerLength > 0)
        {
            node->outer = buildvp(outerPoints, outerLength, d);
        }
    }


    if(innerLength < 1)
    {
        node->inner = NULL;
    }
    if(outerLength < 1)
    {
        node->outer= NULL;
    }


    node->md = median;
    node->vp = point;
    node->idx = id;

    // De-allocate unused memory
    free(points);
    free(distances);
    free(distancesCopy);
    free(innerPoints);
    free(outerPoints);
    free(ids);

    return node;
}

// Return vantage-point subtree with points inside radius
vptree * getInner(vptree * T)
{
    return T->inner;
}

// Return vantage-point subtree with points outside radius
vptree * getOuter(vptree * T)
{
    return T->outer;
}

// Return median of distances to vantage point
double getMD(vptree * T)
{
    return T->md;
}

// Return the coordinates of the vantage point
double * getVP(vptree * T)
{
    return T->vp;
}

// Return the index of the vantage point
int getIDX(vptree * T)
{
    return T->idx;
}

// Calculates the distances of all points from point and writes them to an array. If possible use work-sharing to parallelize
void euclidean(double *point, double *points, double *distances, int n, int d)
{

    // Since
    double accumulator = 0;
    int i;
    if((n-1 > POINT_THRESHOLD) && (PARALLELDIS == true) && (THREADS <= THREADS_MAX - threadCount))
    {
        modThreadCount(THREADS);
        // Note that by removing the numthreads() clause, we can use OpenMP's default value which is usually = #CPU cores
        #pragma omp parallel shared(point, points, distances, n, d) private(i, accumulator) num_threads(THREADS)
        {
            #pragma omp for schedule(static) nowait
            for (i = 0; i < n; i++)
            {
                accumulator = 0;
                for (int j = 0; j < d; j++)
                {
                    accumulator += pow((point[j] - *(points + i * d + j)), 2);
                }
                distances[i] = sqrt(accumulator);
            }
        }
        modThreadCount(-THREADS);
    }else{
        for (i = 0; i < n; i++)
        {
            accumulator = 0;
            for (int j = 0; j < d; j++)
            {
                accumulator += pow((point[j] - *(points + i * d + j)), 2);
            }
            distances[i] = sqrt(accumulator);
        }
    }
    return;
}

// A utility function to swap two elements
void swap(double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

// QuickSort Partition function
// low and high are the range of indexes in arr where partition should work
int partition (double arr[], int low, int high)
{

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
double quickselect_median(double arr[], int length)
{

    if (length % 2 == 1)
    {
        return quickselect(arr, length, (length+1)/2);
    }
    else
    {
        return 0.5 * (quickselect(arr, length, length/2) + quickselect(arr, length, length/2 + 1));
    }
}

// Returns the idx-th element of arr when arr is sorted
// idx is the index (starting from 1) of the point we want to find when the array is sorted. For the median idx should be the middle one (i.e (length+1)/2 for odd lengths etc)
double quickselect(double arr[], int length, int idx)
{

    // Check to end recursion
    if (length == 1)
    {
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
    if (idx <= lowerLength)
    {
        result = quickselect(lower, lowerLength, idx);
    }
    // This means that the median is our pivot point
    else if(idx == pivotIndex)
    {
        result = pivot;
    }
    // This means that the median is in the higher partition
    else
    {
        result = quickselect(higher, higherLength, idx - pivotIndex);
    }

    // Free memory allocated to lower and higher
    free(lower);
    free(higher);

    // Return result
    return result;
}
