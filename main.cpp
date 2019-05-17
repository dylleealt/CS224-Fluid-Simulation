#include "Renderer.h"

int main(){
    Renderer renderer = Renderer(1.0, 60, 500, 100.0, 25.0, 0.0005, 0.2, 0.1, 0.1, 400);
    renderer.simulateAndRender(100,100);

    return 0;
}
