#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define AREA_WIDTH 80
#define AREA_HEIGHT 20
#define NUM_FISH 5
#define MAX_SPEED 2

typedef struct {
    int x;
    int y;
    int speedX;
    int speedY;
} Fish;

void clearScreen() {
    printf("\033[H\033[J");
}

void moveFish(Fish *fish) {
    // Update fish position based on speed
    fish->x += fish->speedX;
    fish->y += fish->speedY;

    // Ensure the fish stays within the defined area
    if (fish->x < 0 || fish->x >= AREA_WIDTH) {
        fish->speedX = -fish->speedX;  // Reverse speed on hitting boundaries
        fish->x += 2 * fish->speedX;   // Move back
    }

    if (fish->y < 0 || fish->y >= AREA_HEIGHT) {
        fish->speedY = -fish->speedY;  // Reverse speed on hitting boundaries
        fish->y += 2 * fish->speedY;   // Move back
    }
}

void drawFish(Fish *fish) {
    printf("><(((('>");  // Fish representation
}

void draw(Fish *fishes, int numFish) {
    clearScreen();

    for (int i = 0; i < numFish; i++) {
        for (int j = 0; j < fishes[i].y; j++) {
            printf("\n");
        }

        for (int j = 0; j < fishes[i].x; j++) {
            printf(" ");
        }

        drawFish(&fishes[i]);
        printf("\n");
    }
}

int main() {
    srand(time(NULL));

    Fish fishes[NUM_FISH];

    // Initialize fish positions and speeds within the defined area
    for (int i = 0; i < NUM_FISH; i++) {
        fishes[i].x = rand() % AREA_WIDTH;
        fishes[i].y = rand() % AREA_HEIGHT;
        fishes[i].speedX = rand() % (2 * MAX_SPEED + 1) - MAX_SPEED;
        fishes[i].speedY = rand() % (2 * MAX_SPEED + 1) - MAX_SPEED;
    }

    while (1) {
        for (int i = 0; i < NUM_FISH; i++) {
            moveFish(&fishes[i]);
        }

        draw(fishes, NUM_FISH);
        usleep(100000);  // Sleep for 100 milliseconds (adjust for desired frame rate)
    }

    return 0;
}
