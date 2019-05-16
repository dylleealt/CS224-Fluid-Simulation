#include "particles.h"

#include <math.h>
#include <iostream>
#include "gl/shaders/ShaderAttribLocations.h"

const glm::vec3 triX = glm::vec3(0.5, 0, 0);
const glm::vec3 triY = glm::vec3(0, 0.5, 0);
const glm::vec3 triZ = glm::vec3(0, 0, 0.5);

Particles::Particles() : m_numRows(50), m_numCols(m_numRows), m_numLayers(m_numRows)
{
}

/**
 * Initializes the Particles by storing positions and normals in a vertex buffer.
 */
void Particles::init() {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Initialize solver
    m_solver = std::make_unique<FluidSolver>();
    m_solver->init(m_numCols, m_numRows, m_numLayers, m_numCols, m_numRows, m_numLayers, 0.0005, 0.2, 0, 0.0167);
    m_solver->reset();
    float visc = m_solver->getViscosity();
    float diff = m_solver->getDiffusionRate();
    float rate = m_solver->getDissipationRate();
    float dt = m_solver->getTimeStep();
    for (int i = 0; i < 1000; ++i){
        std::cout<<i<<std::endl;
        m_solver->update(visc, diff, rate, dt, 3);
    }

    // Initializes a grid of particles, will convert to triangle strip later maybe
    int numParticles = m_numCols * m_numRows * m_numLayers;
    m_particles.reserve(3 * numParticles);
    int index = 0;
    int x_radius = m_numCols / 2;
    int y_radius = m_numRows / 2;
    int z_radius = m_numLayers / 2;
    for (int i = -x_radius; i < x_radius; i+=2){
        for (int j = y_radius; j < m_numRows; j+=2){
            for (int k = -z_radius; k < z_radius; k+=2){
                glm::vec3 pos = glm::vec3(i, j, k);
                m_particles.push_back(pos);
                m_particles.push_back(pos + triY);
                m_particles.push_back(pos + triX);
            }
        }
    }


    // Initialize OpenGLShape.
    m_shape = std::make_unique<OpenGLShape>();
    m_shape->setVertexData(&m_particles[0][0], 3 * m_particles.size(), VBO::GEOMETRY_LAYOUT::LAYOUT_TRIANGLES, m_particles.size());
    m_shape->setAttribute(ShaderAttrib::POSITION, 3, 0, VBOAttribMarker::DATA_TYPE::FLOAT, false);
    m_shape->buildVAO();
}


/**
 * Draws the Particles.
 */
void Particles::draw()
{
//    float visc = m_solver->getViscosity();
//    float diff = m_solver->getDiffusionRate();
//    float rate = m_solver->getDissipationRate();
//    float dt = m_solver->getTimeStep();
//    float hx = m_solver->getHx();
//    float hy = m_solver->getHy();
//    float hz = m_solver->getHz();
//    float **velocity = m_solver->getVelocity();
//    float *vx = velocity[0], *vy = velocity[1], *vz = velocity[2];
//    m_solver->update(visc, diff, rate, dt, 4);
//    // update positions
//    for (int i = 0; i < m_numCols; ++i){
//        for (int j = 0; j < m_numRows; ++j){
//            for (int k = 0; k < m_numLayers; ++k){
//                int index = m_solver->idx(i, j, k);
//                float x = (i + 0.5) * hx;
//                float y = (j + 0.5) * hy;
//                float z = (k + 0.5) * hz;
//                glm::vec3 v = glm::vec3(
//                                m_solver->interpolate(vx, x, y, z),
//                                m_solver->interpolate(vy, x, y, z),
//                                m_solver->interpolate(vz, x, y, z)
//                            );
//                m_particles[index] += v * dt;
//                glm::clamp(m_particles[index], glm::vec3(0), glm::vec3(m_numRows));
//            }
//        }
//    }
    m_shape->draw();
}
