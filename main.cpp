#include "Renderer.h"

int main(){
    Renderer renderer = Renderer(0.0167, 60, 100, 100.0, 25.0, 0.0005, 0.2, 0.1, 0.1, 400);
    renderer.simulateAndRender(300, 300);

    return 0;
}
