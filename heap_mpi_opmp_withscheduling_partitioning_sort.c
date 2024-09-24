#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>  // Include OpenMP header
#include <mpi.h>  // Include MPI header

#define NUM_FISH 100   // Number of fish in the school
#define NUM_STEPS 1000 // Number of simulation steps
#define NUM_THREADS 4   // Number of threads to use (adjust as needed)

#define INITIAL_WEIGHT 1.0            // Initial weight of each fish
#define MAX_WEIGHT_MULTIPLIER 2.0    // Maximum weight multiplier
#define EAT_PROBABILITY 0.2          // Probability of a fish eating (adjust as needed)
#define SWIM_PROBABILITY 0.5         // Probability of a fish swimming (adjust as needed)

// Fish structure
typedef struct {
    double x;
    double y;
    double weight;
} Fish;

// Function to calculate the objective function
double calculateObjectiveFunction(Fish *school, int numFish) {
    double objectiveFunction = 0.0;
    for (int i = 0; i < numFish; i++) {
        double x = school[i].x;
        double y = school[i].y;
        objectiveFunction += sqrt(x * x + y * y);
    }
    return objectiveFunction;
}

// Function to update fish weights based on objective function change
void updateFishWeights(Fish *school, int numFish, double maxDelta) {
    for (int i = 0; i < numFish; i++) {
        school[i].weight += maxDelta;
        if (school[i].weight > INITIAL_WEIGHT * MAX_WEIGHT_MULTIPLIER) {
            school[i].weight = INITIAL_WEIGHT * MAX_WEIGHT_MULTIPLIER;
        }
    }
}

// Function to simulate fish eating (randomly increasing weight)
void eat(Fish *fish) {
    if ((double)rand() / RAND_MAX < EAT_PROBABILITY) {
        double weightIncrease = ((double)rand() / RAND_MAX) * INITIAL_WEIGHT; // Random weight increase
        fish->weight += weightIncrease;
        if (fish->weight > INITIAL_WEIGHT * MAX_WEIGHT_MULTIPLIER) {
            fish->weight = INITIAL_WEIGHT * MAX_WEIGHT_MULTIPLIER;
        }
    }
}

// Function to simulate fish swimming (randomly moving)
void swim(Fish *fish) {
    if ((double)rand() / RAND_MAX < SWIM_PROBABILITY) {
        double deltaX = ((double)rand() / RAND_MAX - 0.5) * 0.1;
        double deltaY = ((double)rand() / RAND_MAX - 0.5) * 0.1;
        fish->x += deltaX;
        fish->y += deltaY;
    }
}

// Function to calculate the barycenter of the fish school
void calculateBarycenter(Fish *school, int numFish, double *barycenterX, double *barycenterY) {
    double sum_qxi_squared_plus_yi_squared = 0.0;

    for (int i = 0; i < numFish; i++) {
        double weight = school[i].weight;
        double qxi_squared_plus_yi_squared = sqrt(school[i].x * school[i].x + school[i].y * school[i].y);

        sum_qxi_squared_plus_yi_squared += qxi_squared_plus_yi_squared;

        *barycenterX += qxi_squared_plus_yi_squared * weight * school[i].x;
        *barycenterY += qxi_squared_plus_yi_squared * weight * school[i].y;
    }

    // Calculate the barycenter coordinates
    if (sum_qxi_squared_plus_yi_squared != 0) {
        *barycenterX /= sum_qxi_squared_plus_yi_squared;
        *barycenterY /= sum_qxi_squared_plus_yi_squared;
    }
}

// Comparison function for sorting fish based on their coordinates
int compareFish(const void *a, const void *b) {
    Fish *fishA = (Fish *)a;
    Fish *fishB = (Fish *)b;

    // Sort based on the x-coordinate, then y-coordinate
    if (fishA->x < fishB->x)
        return -1;
    else if (fishA->x > fishB->x)
        return 1;
    else {
        if (fishA->y < fishB->y)
            return -1;
        else if (fishA->y > fishB->y)
            return 1;
        else
            return 0;
    }
}

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Seed the random number generator based on the MPI rank
    srand(time(NULL) + world_rank);

    // Allocate memory for an array of fish
    int numFishPerProcess = NUM_FISH / world_size;
    Fish *school = (Fish *)malloc(numFishPerProcess * sizeof(Fish));
    if (school == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // Initialize fish randomly within the square
    for (int i = 0; i < numFishPerProcess; i++) {
        school[i].x = ((double)rand() / RAND_MAX - 0.5) * 200;
        school[i].y = ((double)rand() / RAND_MAX - 0.5) * 200;
        school[i].weight = INITIAL_WEIGHT; // Initial weight
    }

    // Sort fish based on coordinates
    qsort(school, numFishPerProcess, sizeof(Fish), compareFish);

    // Record the start time
    double start_time = MPI_Wtime();

    // Perform the simulation for a specified number of steps
    #pragma omp parallel num_threads(NUM_THREADS) // Start OpenMP parallel region
    {
        int numFishPerThread = numFishPerProcess / NUM_THREADS;
        int threadID = omp_get_thread_num();
        int startIdx = threadID * numFishPerThread;
        int endIdx = (threadID == NUM_THREADS - 1) ? numFishPerProcess : startIdx + numFishPerThread;

        // Handle the remainder fish for the last thread
        if (threadID == NUM_THREADS - 1) {
            endIdx = numFishPerProcess;
        }

        // Record the start time for this scheduling
        double scheduling_start_time = MPI_Wtime();

        // Simulate fish for each step
        for (int step = 0; step < NUM_STEPS; step++) {
            // Simulate fish for the current step
            #pragma omp for schedule(static)
            for (int i = startIdx; i < endIdx; i++) {
                eat(&school[i]); // Simulate fish eating
                swim(&school[i]); // Simulate fish swimming
            }

            // Synchronize threads
            #pragma omp barrier

            // Reset the weights for the next simulation step
            #pragma omp for
            for (int i = startIdx; i < endIdx; i++) {
                school[i].weight = INITIAL_WEIGHT;
            }

            // Synchronize threads
            #pragma omp barrier
        }

        // Record the end time for this scheduling and print the execution time for this scheduling
        double scheduling_end_time = MPI_Wtime();
        double scheduling_execution_time = scheduling_end_time - scheduling_start_time;
        printf("Execution Time for Thread %d: %.6lf seconds\n", threadID, scheduling_execution_time);
    } // End of OpenMP parallel region

    // Gather the fish data from all processes
    Fish *allFish = NULL;
    if (world_rank == 0) {
        allFish = (Fish *)malloc(NUM_FISH * sizeof(Fish));
        if (allFish == NULL) {
            printf("Memory allocation failed.\n");
            return 1;
        }
    }

    // Gather fish data using MPI
    MPI_Gather(school, numFishPerProcess, MPI_DOUBLE, allFish, numFishPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculate and print the total execution time in seconds
    double end_time = MPI_Wtime();
    double execution_time = end_time - start_time;
    printf("Total Execution Time: %.6lf seconds\n", execution_time);

    // Free allocated memory for the array of fish
    free(school);

    // Optionally, process and print the final weights if you're rank 0
    if (world_rank == 0) {
        // Process and print the final weights
        printf("Final Weights:\n");
        for (int i = 0; i < NUM_FISH; i++) {
            printf("Fish %d: Weight %.2lf\n", i, allFish[i].weight);
        }

        // Free allocated memory for the array of fish
        free(allFish);
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
