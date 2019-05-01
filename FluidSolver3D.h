#ifndef __FLUIDSOLVER_H__
#define __FLUIDSOLVER_H__

class FluidSolver
{
    public:
        FluidSolver();
        ~FluidSolver();

        void vStep(); // move a vector field forward one time step
        void sStep(); // move a scalar field forward one time step

    private:
      
	void addForce(float *f, float dt, int flag);
        void addSource(float *f, int flag);
        void advect(float *f, float *f0, float **v);
        void diffuse(float *f, float *f0, float k, float dt);
        void project(float *f, float *f0, float dt);
        void dissipate(float *f, float *f0, float a, float dt);
}

#endif
