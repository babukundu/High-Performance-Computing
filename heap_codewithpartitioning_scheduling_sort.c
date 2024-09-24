#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h> // Include OpenMP header

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

int main() {
    // Seed the random number generator
    srand(time(NULL));

    // Allocate memory for an array of fish
    Fish *school = (Fish *)malloc(NUM_FISH * sizeof(Fish));
    if (school == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // Initialize fish randomly within the square
    for (int i = 0; i < NUM_FISH; i++) {
        school[i].x = ((double)rand() / RAND_MAX - 0.5) * 200;
        school[i].y = ((double)rand() / RAND_MAX - 0.5) * 200;
        school[i].weight = INITIAL_WEIGHT; // Initial weight
    }

    // Sort fish based on coordinates
    qsort(school, NUM_FISH, sizeof(Fish), compareFish);

    // Record the start time
    clock_t start_time = clock();

    // Perform the simulation for a specified number of steps
    #pragma omp parallel num_threads(NUM_THREADS) // Start parallel region
    {
        int numFishPerThread = NUM_FISH / NUM_THREADS;
        int threadID = omp_get_thread_num();
        int startIdx = threadID * numFishPerThread;
        int endIdx = (threadID == NUM_THREADS - 1) ? NUM_FISH : startIdx + numFishPerThread;

        // Handle the remainder fish for the last thread
        if (threadID == NUM_THREADS - 1) {
            endIdx = NUM_FISH;
        }

        // Record the start time for this scheduling
        clock_t scheduling_start_time = clock();

        // Print results for each chunk using static scheduling
        printf("Results for Chunk %d using static scheduling:\n", threadID);
        #pragma omp for schedule(static)
        for (int i = startIdx; i < endIdx; i++) {
            eat(&school[i]); // Simulate fish eating
            swim(&school[i]); // Simulate fish swimming
            printf("Fish %d: Weight %.2lf\n", i, school[i].weight);
        }

        // Reset the weights for the next simulation step
        #pragma omp for
        for (int i = startIdx; i < endIdx; i++) {
            school[i].weight = INITIAL_WEIGHT;
        }

        // Record the end time for this scheduling and print the execution time for this scheduling
        clock_t scheduling_end_time = clock();
        double scheduling_execution_time = (double)(scheduling_end_time - scheduling_start_time) / CLOCKS_PER_SEC;
        printf("Execution Time for Chunk %d using static scheduling: %.6lf seconds\n", threadID, scheduling_execution_time);

        // Print results for each chunk using dynamic scheduling
        printf("Results for Chunk %d using dynamic scheduling:\n", threadID);
        scheduling_start_time = clock();
        #pragma omp for schedule(dynamic, 10)
        for (int i = startIdx; i < endIdx; i++) {
            eat(&school[i]); // Simulate fish eating
            swim(&school[i]); // Simulate fish swimming
            printf("Fish %d: Weight %.2lf\n", i, school[i].weight);
        }

        // Reset the weights for the next simulation step
        #pragma omp for
        for (int i = startIdx; i < endIdx; i++) {
            school[i].weight = INITIAL_WEIGHT;
        }

        // Record the end time for this scheduling and print the execution time for this scheduling
        scheduling_end_time = clock();
        scheduling_execution_time = (double)(scheduling_end_time - scheduling_start_time) / CLOCKS_PER_SEC;
        printf("Execution Time for Chunk %d using dynamic scheduling: %.6lf seconds\n", threadID, scheduling_execution_time);

        // Print results for each chunk using guided scheduling
        printf("Results for Chunk %d using guided scheduling:\n", threadID);
        scheduling_start_time = clock();
        #pragma omp for schedule(guided)
        for (int i = startIdx; i < endIdx; i++) {
            eat(&school[i]); // Simulate fish eating
            swim(&school[i]); // Simulate fish swimming
            printf("Fish %d: Weight %.2lf\n", i, school[i].weight);
        }

        // Reset the weights for the next simulation step
        #pragma omp for
        for (int i = startIdx; i < endIdx; i++) {
            school[i].weight = INITIAL_WEIGHT;
        }

        // Record the end time for this scheduling and print the execution time for this scheduling
        scheduling_end_time = clock();
        scheduling_execution_time = (double)(scheduling_end_time - scheduling_start_time) / CLOCKS_PER_SEC;
        printf("Execution Time for Chunk %d using guided scheduling: %.6lf seconds\n", threadID, scheduling_execution_time);
    } // End of parallel region

    // Record the end time
    clock_t end_time = clock();

    // Print the final weights
    printf("Final Weights:\n");
    for (int i = 0; i < NUM_FISH; i++) {
        printf("Fish %d: Weight %.2lf\n", i, school[i].weight);
    }

    // Calculate and print the total execution time in seconds
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Total Execution Time: %.6lf seconds\n", execution_time);

    // Free allocated memory for the array of fish
    free(school);

    return 0;
}
