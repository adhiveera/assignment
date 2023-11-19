#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "coordReader.c"
#include <sys/time.h>  // for wallclock timing functions

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile( int *tour, int tourLength, char *filename);

#define INF 999999

// Define a structure to represent a vertex
typedef struct {
    double x, y;  // Coordinates of the vertex
} Vertex;

int t;
// Function to calculate the Euclidean distance between two vertices
double distance(Vertex v1, Vertex v2) {
    return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y));
}

// Function to find the vertex with the minimal distance from the current vertex
int findMinDistanceVertex(int current, int n, int *visited, Vertex *vertices) {
    int minIndex = -1;
    double minDistance = INF;

    for (int i = 0; i < n; i++) {
        if (!visited[i] && i != current) {
            double dist = distance(vertices[current], vertices[i]);
            if (dist < minDistance) {
                minDistance = dist;
                minIndex = i;
            }
        }
    }

    return minIndex;
}

int *findVk(int current, int n, int *visited, Vertex *vertices, int next) {
    int Vk = -1;
    double minDistance = INF;

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            double dist1 = distance(vertices[current], vertices[i]);
            double dist2 = distance(vertices[next], vertices[i]);
            double total = (dist1 + dist2) - distance(vertices[current], vertices[next]);
            if (total < minDistance) {
                minDistance = total;
                Vk = i;
            }
        }
    }

    static int arr[2];
    arr[0] = next;
    arr[1] = Vk; // next becomes current and the new vertex becomes next.
    return arr;
}

// Function to perform the Cheapest Insertion Algorithm
void *cheapestInsertion(int n, Vertex *vertices) {
    int *visited = (int *)malloc(n * sizeof(int));
    int *tour = (int *)malloc((n + 1) * sizeof(int)); // +1 for returning to the starting point
    int current = 0;                                   // Start with the first vertex
    int next;

    for (int i = 0; i < n; ++i) {
        visited[i] = 0;
    }

    visited[current] = 1; // updating that vertex has been visited
    tour[0] = current;

    // STEP #2
    next = findMinDistanceVertex(current, n, visited, vertices);
    visited[next] = 1; // updating that vertex has been visited

    tour[1] = next;

    // STEP #3
    for (int i = 2; i < n; ++i) {
        int *newVerts = findVk(current, n, visited, vertices, next);
        tour[i] = newVerts[1];
        // insert the new vertex
        current = newVerts[0];
        next = newVerts[1];
        visited[next] = 1;
    }

    // Complete the tour by returning to the starting point
    tour[n] = tour[0];

    // Print the tour
    printf("Tour: ");
    for (int i = 0; i <= n; ++i) {
        printf("%d ", tour[i]);
    }
    printf("\n");
    free(visited);
    return tour;
    
}

int main() {
    /* for timing */
    struct timeval wallStart, wallEnd;

    gettimeofday(&wallStart, NULL); // save start time in to variable 'wallStart'

    int n = readNumOfCoords("9_coords.coord"); // Number of vertices
    Vertex vertices[n];
    double **coordinates = readCoords("9_coords.coord", n);
    
    for (int i = 0; i < n; i++) {
        vertices[i].x = coordinates[i][0];
        vertices[i].y = coordinates[i][1];
    }

    int *t = cheapestInsertion(n, vertices);

    gettimeofday(&wallEnd, NULL); // end time
    double wallSecs = (wallEnd.tv_sec - wallStart.tv_sec);           // just integral number of seconds
    double WALLtimeTaken = 1.0E-06 * ((wallSecs*1000000) + (wallEnd.tv_usec - wallStart.tv_usec)); // and now with any microseconds

    printf("Time taken   =  %f seconds  \n",WALLtimeTaken);
    writeTourToFile(t, n,"result.dat");

    free(t); 
    return 0;
}
