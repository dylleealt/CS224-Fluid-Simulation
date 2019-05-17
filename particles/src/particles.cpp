#include "particles.h"

#include <math.h>
#include <iostream>
#include "gl/shaders/ShaderAttribLocations.h"

const glm::vec3 tri1 = glm::vec3(0.5, 0, 0);
const glm::vec3 tri2 = glm::vec3(0.25, 0, 0.433);
const glm::vec3 tri3 = glm::vec3(0.25, 0.4787, 0.144);

Particles::Particles() : m_numRows(50), m_numCols(m_numRows), m_numLayers(m_numRows)
{
}

/**
 * Initializes the Particles by storing positions and normals in a vertex buffer.
 */
void Particles::init() {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Initializes a grid of particles, will convert to triangle strip later maybe
    int numParticles = m_numCols * m_numRows * m_numLayers;
    int numVerticesPerParticle = 12;
    m_particles.reserve(numVerticesPerParticle * numParticles);
    int index = 0;
    for (int i = 0; i < m_numCols; ++i){
        for (int j = m_numRows / 2; j < m_numRows; ++j){
            for (int k = 0; k < m_numLayers; ++k){
                glm::vec3 pos = glm::vec3(i, j, k);
                // face 1
                m_particles.push_back(pos);
                m_particles.push_back(pos + tri1);
                m_particles.push_back(pos + tri3);
                // face 2
                m_particles.push_back(pos);
                m_particles.push_back(pos + tri2);
                m_particles.push_back(pos + tri1);
                // face 3
                m_particles.push_back(pos);
                m_particles.push_back(pos + tri3);
                m_particles.push_back(pos + tri2);
                // face 4
                m_particles.push_back(pos + tri1);
                m_particles.push_back(pos + tri2);
                m_particles.push_back(pos + tri3);
            }
        }
    }

    // Initialize solver
//    m_solver = std::make_unique<FluidSolver>(m_numCols, m_numRows, m_numLayers, m_numCols, m_numRows, m_numLayers, 0.0005, 0.2, 0, 0.1, 0.0167);
//    m_solver->reset();
//    float visc = m_solver->getViscosity();
//    float diff = m_solver->getDiffusionRate();
//    float rate = m_solver->getDissipationRate();
//    float vorticity = m_solver->getVorticity();
//    float dt = m_solver->getTimeStep();
//    float **velocity = m_solver->getVelocity();
//    float *vx = velocity[0], *vy = velocity[1], *vz = velocity[2];
//    int maxIterations = 200;
//    for (int iteration = 0; iteration < maxIterations; ++iteration){
//        std::cout<<iteration<<std::endl;
//        int flag = (iteration < 0.7f * maxIterations) ? 4 : 1;
//        m_solver->update(visc, diff, rate, vorticity, dt, flag);
//        // update positions
//        for (int i = 0; i < m_particles.size(); i+=numVerticesPerParticle){
//            glm::vec3 pos = m_particles[i];
//            // std::cout<<pos.x<<" "<<pos.y<<" "<<pos.z<<std::endl;
//            glm::vec3 v = glm::vec3(m_solver->interpolate(vx, pos.x, pos.y, pos.z),
//                                    m_solver->interpolate(vy, pos.x, pos.y, pos.z),
//                                    m_solver->interpolate(vz, pos.x, pos.y, pos.z));
//            for (int j = 0; j < numVerticesPerParticle; ++j){
//                m_particles[i + j] += v * dt;
//                glm::clamp(m_particles[index], glm::vec3(1.5), glm::vec3(m_numRows - 1.5));
//            }
//        }
//    }

    // Initialize OpenGLShape.
    m_shape = std::make_unique<OpenGLShape>();
    m_shape->setVertexData(&m_particles[0][0], 3 * m_particles.size(), VBO::GEOMETRY_LAYOUT::LAYOUT_TRIANGLES, m_particles.size());
    m_shape->setAttribute(ShaderAttrib::POSITION, 3, 0, VBOAttribMarker::DATA_TYPE::FLOAT, false);
    m_shape->buildVAO();
    std::cout<<"built VAO!"<<std::endl;
}


/**
 * Draws the Particles.
 */
void Particles::draw()
{
    std::cout<<"draw"<<std::endl;
    m_shape->draw();
    std::cout<<"finished!"<<std::endl;
}
