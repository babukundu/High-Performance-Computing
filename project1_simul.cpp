#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define WIDTH 80
#define HEIGHT 20
#define NUM_OBJECTS 10
#define MAX_VELOCITY 3

typedef struct {
    int x;
    int y;
    int velocityX;
    int velocityY;
} Object;

void update(Object *objects, int numObjects) {
    for (int i = 0; i < numObjects; i++) {
        // Update object position based on velocity
        objects[i].x += objects[i].velocityX;
        objects[i].y += objects[i].velocityY;

        // If object hits the boundary, reverse its velocity
        if (objects[i].x < 0 || objects[i].x >= WIDTH)
            objects[i].velocityX = -objects[i].velocityX;
        if (objects[i].y < 0 || objects[i].y >= HEIGHT)
            objects[i].velocityY = -objects[i].velocityY;
    }
}

void clearScreen() {
    system("clear || cls");
}

void draw(Object *objects, int numObjects) {
    char canvas[HEIGHT][WIDTH];

    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            canvas[i][j] = ' ';
        }
    }

    for (int i = 0; i < numObjects; i++) {
        int x = objects[i].x;
        int y = objects[i].y;
        canvas[y][x] = '*';
    }

    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            printf("%c", canvas[i][j]);
        }
        printf("\n");
    }
}

int main() {
    srand(time(NULL));

    Object objects[NUM_OBJECTS];

    // Initialize objects with random positions and velocities
    for (int i = 0; i < NUM_OBJECTS; i++) {
        objects[i].x = rand() % WIDTH;
        objects[i].y = rand() % HEIGHT;
        objects[i].velocityX = rand() % (2 * MAX_VELOCITY + 1) - MAX_VELOCITY;
        objects[i].velocityY = rand() % (2 * MAX_VELOCITY + 1) - MAX_VELOCITY;
    }

    while (1) {
        clearScreen();
        update(objects, NUM_OBJECTS);
        draw(objects, NUM_OBJECTS);
        _sleep(100000); // Sleep for 100 milliseconds (adjust as needed for desired frame rate)
    }

    return 0;
}
