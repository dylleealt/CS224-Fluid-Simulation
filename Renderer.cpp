#include "Renderer.h"

#include <iostream>
#include <thread>
#include <future>

Renderer::Renderer(float timestep, int numFrames, int numCells, float size, float radius, float visc, float diff, float rate, float vorticity, int numParticles):
    m_numCells(numCells),
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

float Renderer::sdBox(glm::vec3 p, glm::vec3 b)
{
    glm::vec3 d = glm::abs(p) - b;
    return glm::min(glm::max(d.x,glm::max(d.y,d.z)),0.f) + glm::length(glm::max(d,0.f));
}

float Renderer::rayMarch(glm::vec3 rayOrigin, glm::vec3 rayDir) {
    float marchDist = 0.001f;
    float boundingDist = 200.f;
    float threshold = 0.001;
    int maxSteps = 1000;

    for (int i = 0; i < maxSteps; ++i) {
        glm::vec3 p = rayOrigin + (marchDist * rayDir);
        float dist = sdBox(p - glm::vec3(m_numCells / 2), glm::vec3(m_numCells / 2));
        if (dist > boundingDist){
            return -1.f;
        }
//        if (dist < threshold) {
//            return marchDist;
//        }
        if (dist > 0.f){
            marchDist += dist * 0.9f;
        }
        else {
            dist = m_lsSolver->interpolate(m_lsSolver->m_phi, p.x, p.y, p.z);
            if (dist < threshold){
                return marchDist;
            }
            if (dist > boundingDist){
                return marchDist;
            }
            marchDist += dist * 0.9f;
        }
    }

    return -1.f;
}

void Renderer::simulateAndRender(int width, int height) {
    //each iteration of for loop is one step forward in time: define what the timestep is and use this when playing back rendered frames
    for (int frameNum = 0; frameNum < m_numFrames; ++frameNum){
        std::cout<<"Frame: "<<frameNum<<std::endl;
        if (frameNum < 100) {
            m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 4);
        } else {
            m_fSolver->update(m_fSolver->getViscosity(), m_fSolver->getDiffusionRate(), m_fSolver->getDissipationRate(), m_fSolver->getVorticity(), m_timestep, 1);
        }
       m_lsSolver->update(m_fSolver->getVelocity());
//       float *phi = m_lsSolver->m_phi;
//       for (int i = 0; i < 100 * 100 * 100; ++i){
//           std::cout<<phi[i]<<std::endl;
//       }

       QImage image(width, height, QImage::Format_RGB32);
       QRgb *imageData = reinterpret_cast<QRgb *>(image.bits());

       glm::vec3 cameraPos = glm::vec3(-50, 60, -50);
       glm::vec3 target = glm::vec3(0.0, 30.0, 0.0);
       glm::vec3 lookVector = glm::normalize(cameraPos - target);
       glm::vec3 up = glm::vec3(0.0, 1.0, 0.0);

       float focalLength = 2.0;
       glm::vec3 camFwrd = lookVector * -1.f;
       glm::vec3 camRght = glm::normalize(glm::cross(camFwrd, up));
       glm::vec3 camUp = glm::normalize(glm::cross(camRght, camFwrd));
       std::cout<<"Rendering..."<<std::endl;
       int num_pixels = width * height;
//       // MULTITHREADING
//       int num_threads = std::thread::hardware_concurrency();
//       // set up future vector
//       std::vector<std::future<void>> results;
//       for (int i = 0; i < num_threads; i++){ results.emplace_back(std::async([=](){
//       for (int offset = i; offset < num_pixels; offset += num_threads){
         for (int offset = 0; offset < num_pixels; ++offset){
               int x = offset % width;
               int y = offset / width;
               glm::vec3 rayDir((2.f * x / width) - 1, 1 - (2.f * y / height), focalLength); // this is in camera coords, not world coords. Can I use the projectOnto function from before to shift the coordinates???
               rayDir = glm::normalize(camRght * rayDir.x + camUp * rayDir.y + camFwrd * rayDir.z);
               float dist = rayMarch(cameraPos, rayDir);
               if (dist < 0.f) {
                   imageData[offset] = qRgb(0.0, 0.0, 0.0);
               } else {
                   int iDist = (int) dist;
//                   std::cout<<iDist<<std::endl;
                   imageData[offset] = qRgb(iDist, iDist, iDist);
               }
           }
       // })); }
       QString output = QString("../img").append(QString::number(frameNum)).append(QString(".png"));
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
