#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include "../inc/vptree.h"

#define THREADS 2

// Function Prototypes
vptree * buildvp(double *X, int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double getMD(vptree * T);
double * getVP(vptree * T);
int getIDX(vptree * T);
void *euclidean(void *arg);
void swap(double *a, double *b);
int partition (double arr[], int low, int high);
double quickselect_median(double arr[], int length);
double quickselect(double arr[], int length, int idx);
void *buildvp_wrapper(void *arg);

// Flag used to detect if build_vp has already been called. If it has not, X is the original input array.
// If it has it means that X is the points vector with an idx vector extended to it at the end
bool runFlag = false;
// Development flags to switch execution mode from serial to parallel for distance calculation and subtree creation
bool parallelDis = true;
bool parallelSub = true;

// Struct used to pass arguments to threads for distance calculation
typedef struct
{
    int tid, n, d;
    double *point, *points, *distances;
} dtargs;

// Struct used to pass arguments to threads for subtree creation
typedef struct
{
    int tid, n, d;
    double *X;
    vptree *subtree;
} stargs;

// Wrapper function to use buildvp() with pthreads
void *buildvp_wrapper(void *arg)
{
    ((stargs *)arg)->subtree = buildvp(((stargs *)arg)->X, ((stargs *)arg)->n, ((stargs *)arg)->d);
    return;
}

// Function that recursively builds the binary tree
vptree * buildvp(double *X, int n, int d)
{

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

    // Block size * threads = points number (n-1)
    pthread_t disThread[THREADS];
    dtargs disArg[THREADS];
    // blockSize = Points per thread
    int blockSize = floor((float)(n-1) / THREADS);

    //TODO: Point number assignment per thread bellow can probably be done in a better way
    // Calculate distances in parallel if true, else do it sequentially
    if(parallelDis == true)
    {
        // Create threads
        for (int i = 0; i < THREADS; i++)
        {
            disArg[i].d = d;
            disArg[i].point = point;
            disArg[i].tid = i;
            disArg[i].points = points + i * blockSize * d;
            disArg[i].distances = distances + i * blockSize;

            // Check how many points to assign to each block. The last block gets the remaining points
            if( i < THREADS - 1)
            {
                disArg[i].n = blockSize;
            }
            else
            {
                disArg[i].n = (n-1) - blockSize * i;
            }
            pthread_create(&disThread[i], NULL, euclidean, (void *)&disArg[i]);
        }

        // Join threads
        for (int i = 0; i < THREADS; i++)
        {
            pthread_join(disThread[i], NULL);
        }
    }
    else
    {
        disArg[0].d = d;
        disArg[0].point = point;
        disArg[0].tid = 0;
        disArg[0].points = points;
        disArg[0].n = n-1;
        disArg[0].distances = distances;
        euclidean((void *)&disArg[0]);
    }


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

    // Booleans to keep track whether a thread has been created to work on a subtree
    bool threadActive[2];

    // Thread and stargs arrays for parallel subtree threads
    pthread_t subThread[2];
    stargs subArg[2];

    // Build subtrees in parallel or sequentially
    if(parallelSub == true)
    {

        // Create threads
        if(innerLength > 0)
        {
            threadActive[0] = true;
            subArg[0].d = d;
            subArg[0].n = innerLength;
            subArg[0].subtree = node->inner;
            subArg[0].tid = 0;
            subArg[0].X = innerPoints;
            pthread_create(&subThread[0], NULL, buildvp_wrapper, (void *)&subArg[0]);
        }

        if(outerLength > 0)
        {
            threadActive[1] = true;
            subArg[1].d = d;
            subArg[1].n = outerLength;
            subArg[1].subtree = node->outer;
            subArg[1].tid = 1;
            subArg[1].X = outerPoints;
            pthread_create(&subThread[1], NULL, buildvp_wrapper, (void *)&subArg[1]);
        }

        // Join threads
        for(int i=0; i<2; i++)
        {
            if(threadActive[i] == true)
            {
                pthread_join(subThread[i], NULL);
            }
        }

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

// Calculates the distances of all points from point and writes them to an array
void *euclidean(void *arg)
{

    // Retrieve variables by casting the argument to a dtarg struct pointer
    double *point = ((dtargs *) arg) -> point;
    double *points = ((dtargs *) arg) -> points;
    double *distances = ((dtargs *) arg) -> distances;
    int n = ((dtargs *) arg) -> n;
    int d = ((dtargs *) arg) -> d;

    double accumulator = 0;

    for (int i = 0; i < n; i++)
    {
        accumulator = 0;
        for (int j = 0; j < d; j++)
        {
            accumulator += pow((point[j] - *(points + i * d + j)), 2);
        }
        distances[i] = sqrt(accumulator);
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

