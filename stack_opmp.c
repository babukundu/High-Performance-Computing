#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>  // Added OpenMP header

#define NUM_FISH 100   // Number of fish in the school
#define NUM_STEPS 2500 // Number of simulation steps
#define NUM_THREADS 4   // Number of threads to use (modify as needed)
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

// Fish action types
typedef enum {
    EAT_ACTION,
    SWIM_ACTION
} FishActionType;

// Fish action
typedef struct {
    FishActionType actionType;
    int fishIndex;
} FishAction;

// Stack structure for fish actions
typedef struct {
    FishAction *actions;
    int top;
    int capacity;
} ActionStack;

// Function to initialize the action stack
ActionStack *initActionStack(int capacity) {
    ActionStack *stack = (ActionStack *)malloc(sizeof(ActionStack));
    if (stack == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    stack->actions = (FishAction *)malloc(capacity * sizeof(FishAction));
    if (stack->actions == NULL) {
        free(stack);
        printf("Memory allocation failed.\n");
        exit(1);
    }

    stack->top = -1;
    stack->capacity = capacity;

    return stack;
}

// Function to push an action onto the stack
void pushAction(ActionStack *stack, FishAction action) {
    if (stack->top == stack->capacity - 1) {
        printf("Stack overflow.\n");
        exit(1);
    }

    stack->top++;
    stack->actions[stack->top] = action;
}

// Function to pop an action from the stack
FishAction popAction(ActionStack *stack) {
    if (stack->top == -1) {
        printf("Stack underflow.\n");
        exit(1);
    }

    FishAction action = stack->actions[stack->top];
    stack->top--;
    return action;
}

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

    // Create an action stack to track fish actions
    ActionStack *actionStack = initActionStack(NUM_FISH * NUM_STEPS);

    // Record the start time
    clock_t start_time = clock();

    // Set the number of threads
    omp_set_num_threads(NUM_THREADS);

    // Perform the simulation for a specified number of steps
    for (int step = 0; step < NUM_STEPS; step++) {
#pragma omp parallel for
        for (int i = 0; i < NUM_FISH; i++) {
            FishAction fishAction;
            fishAction.fishIndex = i;

            if ((double)rand() / RAND_MAX < EAT_PROBABILITY) {
                fishAction.actionType = EAT_ACTION;
            } else {
                fishAction.actionType = SWIM_ACTION;
            }

            pushAction(actionStack, fishAction);
        }

        // Process fish actions from the stack in parallel
#pragma omp parallel for
        for (int i = 0; i <= actionStack->top; i++) {
            FishAction currentAction = actionStack->actions[i];
            int fishIndex = currentAction.fishIndex;

            if (currentAction.actionType == EAT_ACTION) {
                eat(&school[fishIndex]);
            } else if (currentAction.actionType == SWIM_ACTION) {
                swim(&school[fishIndex]);
            }
        }

        // Calculate the barycenter in parallel
        double barycenterX = 0.0;
        double barycenterY = 0.0;

#pragma omp parallel for reduction(+:barycenterX, barycenterY)
        for (int i = 0; i < NUM_FISH; i++) {
            double weight = school[i].weight;
            double qxi_squared_plus_yi_squared = sqrt(school[i].x * school[i].x + school[i].y * school[i].y);

            barycenterX += qxi_squared_plus_yi_squared * weight * school[i].x;
            barycenterY += qxi_squared_plus_yi_squared * weight * school[i].y;
        }

        // Update fish weights based on the change in the objective function in parallel
        double maxDelta = 0.0;
#pragma omp parallel for reduction(max:maxDelta)
        for (int i = 0; i < NUM_FISH; i++) {
            double originalX = school[i].x;
            double originalY = school[i].y;

            // Simulate fish swimming
            double deltaX = ((double)rand() / RAND_MAX - 0.5) * 0.2;
            double deltaY = ((double)rand() / RAND_MAX - 0.5) * 0.2;
            school[i].x += deltaX;
            school[i].y += deltaY;

            // Calculate change in objective function
            double deltaObjective = calculateObjectiveFunction(school, NUM_FISH) - calculateObjectiveFunction(school, NUM_FISH);

            // Update maximum delta
            if (deltaObjective > maxDelta) {
                maxDelta = deltaObjective;
            }

            // Restore fish position
            school[i].x = originalX;
            school[i].y = originalY;
        }

        // Update fish weights
#pragma omp parallel for
        for (int i = 0; i < NUM_FISH; i++) {
            school[i].weight += maxDelta;
            if (school[i].weight > INITIAL_WEIGHT * MAX_WEIGHT_MULTIPLIER) {
                school[i].weight = INITIAL_WEIGHT * MAX_WEIGHT_MULTIPLIER;
            }
        }
    }

    // Record the end time
    clock_t end_time = clock();

    // Calculate and print the final objective function
    double finalObjective = calculateObjectiveFunction(school, NUM_FISH);
    printf("Final Objective Function: %lf\n", finalObjective);

    // Free allocated memory for the array of fish and the action stack
    free(school);
    free(actionStack->actions);
    free(actionStack);

    // Calculate and print the execution time in seconds
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution Time: %.6lf seconds\n", execution_time);

    return 0;
}
