#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

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
void updateFishWeights(Fish *school, int numFish) {
    double maxDelta = 0.0;
    for (int i = 0; i < numFish; i++) {
        double originalX = school[i].x;
        double originalY = school[i].y;

        // Simulate fish swimming
        double deltaX = ((double)rand() / RAND_MAX - 0.5) * 0.2;
        double deltaY = ((double)rand() / RAND_MAX - 0.5) * 0.2;
        school[i].x += deltaX;
        school[i].y += deltaY;

        // Calculate change in objective function
        double deltaObjective = calculateObjectiveFunction(school, numFish) - calculateObjectiveFunction(school, numFish);

        // Update maximum delta
        if (deltaObjective > maxDelta) {
            maxDelta = deltaObjective;
        }

        // Restore fish position
        school[i].x = originalX;
        school[i].y = originalY;
    }

    // Update fish weights
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

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Seed the random number generator
    srand(time(NULL));

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
    double start_time, end_time;
    double execution_time;

    // Perform the simulation for a specified number of steps

    // Static Scheduling
    if (rank == 0) {
        printf("Static Scheduling:\n");
        start_time = MPI_Wtime();
    }

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for schedule(static)
        for (int step = 0; step < NUM_STEPS; step++) {
            for (int i = 0; i < numFishPerProcess; i++) {
                eat(&school[i]);
                swim(&school[i]);
            }

            double barycenterX = 0.0;
            double barycenterY = 0.0;
            calculateBarycenter(school, numFishPerProcess, &barycenterX, &barycenterY);

            updateFishWeights(school, numFishPerProcess);
        }
    }

    if (rank == 0) {
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;
        printf("Final Objective Function: %lf\n", calculateObjectiveFunction(school, numFishPerProcess));
        printf("Execution Time: %.6lf seconds\n", execution_time);
    }

    // Dynamic Scheduling
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before starting dynamic scheduling

    if (rank == 0) {
        printf("\nDynamic Scheduling:\n");
        start_time = MPI_Wtime();
    }

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for schedule(dynamic)
        for (int step = 0; step < NUM_STEPS; step++) {
            for (int i = 0; i < numFishPerProcess; i++) {
                eat(&school[i]);
                swim(&school[i]);
            }

            double barycenterX = 0.0;
            double barycenterY = 0.0;
            calculateBarycenter(school, numFishPerProcess, &barycenterX, &barycenterY);

            updateFishWeights(school, numFishPerProcess);
        }
    }

    if (rank == 0) {
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;
        printf("Final Objective Function: %lf\n", calculateObjectiveFunction(school, numFishPerProcess));
        printf("Execution Time: %.6lf seconds\n", execution_time);
    }

    // Guided Scheduling
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before starting guided scheduling

    if (rank == 0) {
        printf("\nGuided Scheduling:\n");
        start_time = MPI_Wtime();
    }

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for schedule(guided)
        for (int step = 0; step < NUM_STEPS; step++) {
            for (int i = 0; i < numFishPerProcess; i++) {
                eat(&school[i]);
                swim(&school[i]);
            }

            double barycenterX = 0.0;
            double barycenterY = 0.0;
            calculateBarycenter(school, numFishPerProcess, &barycenterX, &barycenterY);

            updateFishWeights(school, numFishPerProcess);
        }
    }

    if (rank == 0) {
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;
        printf("Final Objective Function: %lf\n", calculateObjectiveFunction(school, numFishPerProcess));
        printf("Execution Time: %.6lf seconds\n", execution_time);
    }

    // Free allocated memory for the array of fish
    free(school);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
