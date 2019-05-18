#ifndef PARTICLES_H
#define PARTICLES_H

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"            // glm::vec*, mat*, and basic glm functions
#include "glm/gtx/transform.hpp"  // glm::translate, scale, rotate
#include "glm/gtc/type_ptr.hpp"   // glm::value_ptr
#include <vector>

#include "openglshape.h"
#include "../FluidSolver.h"

class Particles {
public:
    Particles();

    void init();
    void draw();
    void test();
    int getFrameNumber() { return m_frameNumber; }

private:
    float randValue(int row, int col);
    const float m_numRows, m_numCols, m_numLayers;

    std::unique_ptr<OpenGLShape> m_shape;
    std::unique_ptr<FluidSolver> m_solver;

    int m_frameNumber;

    std::vector<glm::vec3> m_particles;
};

#endif // PARTICLES_H
