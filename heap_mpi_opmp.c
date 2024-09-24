#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define NUM_FISH 100000
#define NUM_THREADS 4

typedef struct {
    double x;
    double y;
    double weight;
} Fish;

void generateFishData(Fish *fish, int numFish) {
    srand(time(NULL));
    for (int i = 0; i < numFish; i++) {
        fish[i].x = ((double)rand() / RAND_MAX - 0.5) * 200;
        fish[i].y = ((double)rand() / RAND_MAX - 0.5) * 200;
        fish[i].weight = 1.0;
    }
}

void writeFishDataToFile(const char *filename, Fish *fish, int numFish) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing: %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < numFish; i++) {
        fprintf(file, "%f %f %f\n", fish[i].x, fish[i].y, fish[i].weight);
    }

    fclose(file);
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time, end_time;

    int fishPerProcess = NUM_FISH / size;
    int remainingFish = NUM_FISH % size;

    Fish *fish = (Fish *)malloc(fishPerProcess * sizeof(Fish));
    if (fish == NULL) {
        printf("Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Record the start time
    start_time = MPI_Wtime();

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        int chunkSize = fishPerProcess / num_threads;
        int chunkStart = tid * chunkSize;
        int chunkEnd = (tid == num_threads - 1) ? fishPerProcess : chunkStart + chunkSize;

        generateFishData(&fish[chunkStart], chunkEnd - chunkStart);
    }

    Fish *recvFish = (Fish *)malloc(fishPerProcess * sizeof(Fish));
    if (recvFish == NULL) {
        printf("Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Distribute fish data from the master to all processes
    MPI_Scatter(fish, fishPerProcess, MPI_DOUBLE, recvFish, fishPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Gather fish data back to the master
    MPI_Gather(recvFish, fishPerProcess, MPI_DOUBLE, fish, fishPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Record the end time
    end_time = MPI_Wtime();

    // Write fish data to a file immediately after the master generates it
    writeFishDataToFile("fish_data_before.txt", fish, fishPerProcess + remainingFish);

    // Write fish data to another file after the master gets back all the data
    if (rank == 0) {
        writeFishDataToFile("fish_data_after.txt", fish, fishPerProcess + remainingFish);
    }

    free(fish);
    free(recvFish);

    // Calculate and print the total execution time
    if (rank == 0) {
        printf("Total execution time: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}