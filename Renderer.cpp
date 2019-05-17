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
                if (dist < half_diag) {
                    m_fSolver->m_d0[x + numCells * (y + numCells * z)] = 1.0;
                } else {
//                    std::cerr<<"test\n";
                    m_fSolver->m_d0[x + numCells * (y + numCells * z)] = 0.0;
                }
            }
        }
    }
}

Renderer::~Renderer()
{
}

////TODO: combine this with density field for efficiency! Also combine lsSolver and fSolver
//static float* initializeDistanceField(int numCells, float size, float radius) {
//    float cellSize = size/numCells;
//    glm::vec3 center = glm::vec3(size*0.5,size*0.5,size*0.5);
//    float* phi = new float[numCells*numCells*numCells];
//    for (int x = 0; x < numCells; x++) {
//        for (int y = 0; y < numCells; y++) {
//            for (int z = 0; z < numCells; z++) {
//                float dist = glm::distance(cellSize * glm::vec3(x, y, z), center) - radius;
//                phi[x + numCells * (y + numCells * z)] = dist;
//            }
//        }
//    }
//    return phi;
//}

//static float* initializeDensityField(int numCells, float size, float radius) {
//    float cellSize = size/numCells;
//    float half_diag = sqrt(cellSize*cellSize*3.0)*0.5;
//    glm::vec3 center = glm::vec3(size*0.5,size*0.5,size*0.5);
//    float* density = new float[numCells*numCells*numCells];
//    for (int x = 0; x < numCells; x++) {
//        for (int y = 0; y < numCells; y++) {
//            for (int z = 0; z < numCells; z++) {
//                float dist = glm::distance(cellSize * glm::vec3(x, y, z), center) - radius;
//                if (dist < half_diag) {
//                    density[x + numCells * (y + numCells * z)] = 1.0;
//                } else {
//                    density[x + numCells * (y + numCells * z)] = 0.0;
//                }
//            }
//        }
//    }
//    return density;
//}

float Renderer::rayMarch(glm::vec3 rayOrigin, glm::vec3 rayDir) {
    float distance_to_origin = 0.0;
    int max_steps = 200;
    float max_dist = 100;
    float surface_dist = 0.005;
    for (int i = 0; i < max_steps; i++) {
        glm::vec3 p = rayOrigin + (distance_to_origin * rayDir);
        float distance_to_surface = m_lsSolver->interpolate(m_lsSolver->m_phi, p.x, p.y, p.z);
        distance_to_origin += distance_to_surface;
        if(distance_to_surface < surface_dist || distance_to_origin > max_dist) {
            break;
        }
    }
    return distance_to_origin;
}

glm::mat3x3 Renderer::relativeMatrix(glm::vec3 lookVector) {
    float z = lookVector.z;
    if (z == 1.0) {
        return glm::mat3x3(1.0, 0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    }
    if (z == -1.0) {
        glm::mat3x3 toReturn(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
//        toReturn << (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
        return toReturn;
    }
    float x = lookVector.x;
    float y = lookVector.y;

    float xx_yy = (x * x) + (y * y);
    assert (xx_yy > 0.0);
    glm::mat3x3 toReturn(((x * x * z) + (y * y)) / xx_yy,  (y * y) + ((x * y * z) / xx_yy),  x,
                          (x * y * (z - 1.0)) / xx_yy,     ((x * x) + (y * y * z)) / xx_yy,  y,
                          -x,                              -y,                               z);
    return toReturn;
}

void Renderer::simulateAndRender(int width, int height) {

    float** frames;
    //each iteration of for loop is one step forward in time: define what the timestep is and use this when playing back rendered frames
    for (int i = 0; i < m_numFrames; ++i){
        if (i < 100) {
            m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 4);
        } else {
            m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 1);
        }
       m_lsSolver->update(m_fSolver->getVelocity());
       QImage image(width, height, QImage::Format_RGB32);
       QRgb *imageData = reinterpret_cast<QRgb *>(image.bits());
       glm::vec3 cameraPos = glm::vec3(0,1,0);
       glm::vec3 lookVector = glm::vec3(0,0,0) - cameraPos; //TODO: triple check this because it's definitely going to be wrong
       for(int y = 0; y < height; ++y) {
           for(int x = 0; x < width; ++x) {
               int offset = x + (y * width);
               glm::vec3 rayDir((2.f * x / width) - 1, 1 - (2.f * y / height), -1); // this is in camera coords, not world coords. Can I use the projectOnto function from before to shift the coordinates???
               glm::normalize(rayDir);
               rayDir = rayDir * relativeMatrix(lookVector);
               float dist = rayMarch(cameraPos, rayDir);
               int iDist = (int) dist;
               imageData[offset] = qRgb(iDist, iDist, iDist);
           }
       }
       QString output = QString("C:\\Users\\Thomas\\Documents\\Brown\\cs2240\\images\\img").append(QString::number(i)).append(QString(".png"));
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
