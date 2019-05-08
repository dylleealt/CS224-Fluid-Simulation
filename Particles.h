#ifndef __PARTICLES_H__
#define __PARTICLES_H__

inline int idx(int i, int j, int k){ return i + width * (j + height * k); }

class Particles
{
    public:
        Particles();
        ~Particles();
        void reset();

        int getNumParticles(){ return numParticles; }
       
        float getTimeStep(){ return dt; }

		glm::vec3 *getPositions(){ return positions; }
		glm::vec3 interpolate_velocity(glm::vec3 pos);
		void updatePositions();

    private:
		int numParticles; // = 1000000; // number of particles
		float dt; // THIS SHOULD BE THE SAME AS THE FLUID GRID TIMESTEP!!
		glm::vec3 positions[];
}

#endif
