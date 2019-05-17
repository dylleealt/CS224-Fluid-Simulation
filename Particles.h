#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include <memory>
#include <vector>

#include "glm.hpp"

class Particles
{
    public:
        Particles();
        ~Particles();

//        void initialize_particles();

//        int getNumParticles(){ return numParticles; }
       
//        float getTimeStep(){ return dt; }

//		glm::vec3 *getPositions(){ return positions; }
//        glm::vec3 interpolate_velocity(glm::vec3 pos, float**v);
//		void updatePositions();

    private:
//		int numParticles; // = 1000000; // number of particles
//		float dt; // THIS SHOULD BE THE SAME AS THE FLUID GRID TIMESTEP!!
//		glm::vec3 positions[];
};

#endif
