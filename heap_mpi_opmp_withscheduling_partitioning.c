#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h> // Include OpenMP header
#include <mpi.h> // Include MPI header

#define NUM_FISH 100   // Number of fish in the school
#define NUM_STEPS 1000 // Number of simulation steps
#define NUM_THREADS 4   // Number of threads to use (adjust as needed)

#define INITIAL_WEIGHT 1.0 // Initial weight of each fish
#define MAX_WEIGHT_MULTIPLIER 2.0 // Maximum weight multiplier
#define EAT_PROBABILITY 0.2 // Probability of a fish eating (adjust as needed)
#define SWIM_PROBABILITY 0.5 // Probability of a fish swimming (adjust as needed)

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

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Seed the random number generator
    srand(time(NULL) + rank); // Use different seeds for each MPI process

    // Determine the number of fish to simulate in each MPI process
    int numFishPerProcess = NUM_FISH / size;

    // Allocate memory for an array of fish
    Fish *school = (Fish *)malloc(numFishPerProcess * sizeof(Fish));
    if (school == NULL) {
        printf("Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize fish randomly within the square
    for (int i = 0; i < numFishPerProcess; i++) {
        school[i].x = ((double)rand() / RAND_MAX - 0.5) * 200;
        school[i].y = ((double)rand() / RAND_MAX - 0.5) * 200;
        school[i].weight = INITIAL_WEIGHT; // Initial weight
    }

    // Record the start time
    double start_time = MPI_Wtime();

    // Perform the simulation for a specified number of steps

    // Parallel region using OpenMP threads within each MPI process
    #pragma omp parallel num_threads(NUM_THREADS)
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
        double scheduling_start_time = omp_get_wtime();

        // Simulation steps
        for (int step = 0; step < NUM_STEPS; step++) {
            // Fish behavior simulation
            for (int i = startIdx; i < endIdx; i++) {
                eat(&school[i]);
                swim(&school[i]);
            }

            // Synchronize threads to ensure data consistency
            #pragma omp barrier

            // Calculate and share maxDelta among threads
            double maxDelta = 0.0;
            for (int i = startIdx; i < endIdx; i++) {
                double originalX = school[i].x;
                double originalY = school[i].y;

                // Calculate change in objective function
                double deltaObjective = calculateObjectiveFunction(school, numFishPerProcess) - calculateObjectiveFunction(school, numFishPerProcess);

                // Update maximum delta
                if (deltaObjective > maxDelta) {
                    maxDelta = deltaObjective;
                }

                // Restore fish position
                school[i].x = originalX;
                school[i].y = originalY;
            }

            // Synchronize threads to ensure maxDelta calculation is complete
            #pragma omp barrier

            // Update fish weights based on maxDelta
            updateFishWeights(school, numFishPerProcess, maxDelta);

            // Synchronize threads to ensure weights are updated
            #pragma omp barrier
        }

        // Record the end time for this scheduling and print the execution time for this scheduling
        double scheduling_end_time = omp_get_wtime();
        printf("Thread %d Execution Time: %.6lf seconds\n", threadID, scheduling_end_time - scheduling_start_time);
    }

    // Synchronize MPI processes
    MPI_Barrier(MPI_COMM_WORLD);

    // Gather fish data from all processes into the root process (rank 0)
    Fish *allSchool = NULL;
    if (rank == 0) {
        allSchool = (Fish *)malloc(NUM_FISH * sizeof(Fish));
        if (allSchool == NULL) {
            printf("Memory allocation failed in the root process.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Gather(school, numFishPerProcess * sizeof(Fish), MPI_BYTE, allSchool, numFishPerProcess * sizeof(Fish), MPI_BYTE, 0, MPI_COMM_WORLD);

    // The root process (rank 0) has all the data
    if (rank == 0) {
        // Record the end time
        double end_time = MPI_Wtime();

        // Print the final objective function and execution time
        double finalObjectiveFunction = calculateObjectiveFunction(allSchool, NUM_FISH);
        printf("Final Objective Function: %lf\n", finalObjectiveFunction);
        printf("Total Execution Time: %.6lf seconds\n", end_time - start_time);

        // Free allocated memory for the array of fish
        free(allSchool);
    }

    // Free allocated memory for the array of fish in each MPI process
    free(school);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
