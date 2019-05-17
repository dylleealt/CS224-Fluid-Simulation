#include "Renderer.h"

#include <iostream>

Renderer::Renderer(float timestep, int numFrames, int numCells, float size, float radius, float visc, float diff, float rate, float vorticity, int numParticles):
    m_timestep(timestep),
    m_numFrames(numFrames)
{
    m_fSolver = std::make_unique<FluidSolver>(numCells, numCells, numCells, size, size, size, visc, diff, rate, vorticity, timestep);
    m_lsSolver = std::make_unique<LevelSetSolver>(numCells, numCells, numCells, size, size, size, timestep, numParticles);

    float cellSize = size/numCells;
    float half_diag = sqrt(cellSize*cellSize*3.0)*0.5;
    glm::vec3 center = glm::vec3(size*0.5,size*0.5,size*0.5);
    for (int x = 0; x < numCells; x++) {
        for (int y = 0; y < numCells; y++) {
            for (int z = 0; z < numCells; z++) {
                float dist = glm::distance(cellSize * glm::vec3(x, y, z), center) - radius;
                m_lsSolver->m_phi0[x + numCells * (y + numCells * z)] = dist;
                m_lsSolver->m_phi[x + numCells * (y + numCells * z)] = dist;
                if (dist < half_diag) {
                    m_fSolver->m_d0[x + numCells * (y + numCells * z)] = 1.0;
                } else {
                    m_fSolver->m_d0[x + numCells * (y + numCells * z)] = 0.0;
                }
            }
        }
    }
//    m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 4);
//    m_lsSolver->update(m_fSolver->getVelocity());
}

Renderer::~Renderer()
{
}

float Renderer::rayMarch(glm::vec3 rayOrigin, glm::vec3 rayDir) {
    float distance_to_origin = 0.0;
    int max_steps = 200;
    float max_dist = 1000;
    float surface_dist = 0.005;
    for (int i = 0; i < max_steps; i++) {
//        std::cerr<<"test1\n";
        glm::vec3 p = rayOrigin + (distance_to_origin * rayDir);
        std::cerr<<"position along ray: "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
        if(p.x >= m_lsSolver->getNx() || p.x <= 1 || p.y >= m_lsSolver->getNy() || p.z <= 1 || p.x >= m_lsSolver->getNz() || p.z <= 1) { //raymarch to cube then use distance field inside
             std::cerr<<"test2\n";
            return 0.0;
            break;
        }
//        std::cerr<<"test3: "<<p.x<<", "<<p.y<<", "<<p.z<<"\n";
        float distance_to_surface = m_lsSolver->interpolate(m_lsSolver->m_phi, p.x, p.y, p.z);
        distance_to_origin += distance_to_surface;
        if(distance_to_surface < surface_dist || distance_to_origin > max_dist) {

            if(distance_to_origin > max_dist) {
                std::cerr<<"test0\n";
//                distance_to_origin = INFINITY;
                return 0.0;
            }
            return 100000 * distance_to_origin;
//            break;
        }
    }
    return distance_to_origin;
}

void Renderer::simulateAndRender(int width, int height) {
//    std::cerr<<"test\n";

//    float** frames;
    //each iteration of for loop is one step forward in time: define what the timestep is and use this when playing back rendered frames
    for (int i = 0; i < 1; ++i){
//    for (int i = 0; i < m_numFrames; ++i){
        if (i < 100) {
            m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 4);
        } else {
            m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 1);
        }
       m_lsSolver->update(m_fSolver->getVelocity());
       QImage image(width, height, QImage::Format_RGB32);
       QRgb *imageData = reinterpret_cast<QRgb *>(image.bits());
       glm::vec3 cameraPos = glm::vec3(98.0,98.0,98.0);
       glm::vec3 target = glm::vec3(0.0,0.0,0.0);
       glm::vec3 lookVector = glm::normalize(cameraPos - target);
       glm::vec3 up = glm::vec3(0.0, 1.0, 0.0);

       glm::vec3 camFwrd = lookVector * -1.f;
       glm::vec3 camRght = glm::normalize(glm::cross(camFwrd, up));
       glm::vec3 camUp = glm::normalize(glm::cross(camRght, camFwrd));
       for(int y = 0; y < height; ++y) {
           for(int x = 0; x < width; ++x) {
               int offset = x + (y * width);
               glm::vec3 rayDir((2.f * x / width) - 1, 1 - (2.f * y / height), -1); // this is in camera coords, not world coords. Can I use the projectOnto function from before to shift the coordinates???
               rayDir = glm::normalize(camRght * rayDir.x + camUp * rayDir.y + camFwrd * rayDir.z);
               float dist = rayMarch(cameraPos, rayDir);
               std::cerr<<dist<<"\n";
               int iDist = (int) dist;
               imageData[offset] = qRgb(x, y, iDist);
           }
       }
       QString output = QString("../img").append(QString::number(i)).append(QString(".png"));
       bool success = image.save(output);
       if(!success) {
           success = image.save(output, "PNG");
       }
       if(success) {
           std::cout << "Wrote rendered image to " << output.toStdString() << std::endl;
       } else {
           std::cerr << "Error: failed to write image to " << output.toStdString() << std::endl;
       }
    }
}
