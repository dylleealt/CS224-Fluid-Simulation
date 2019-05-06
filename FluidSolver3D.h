#ifndef __FLUIDSOLVER_H__
#define __FLUIDSOLVER_H__

class FluidSolver
{
    public:
        FluidSolver();
        ~FluidSolver();

        void vStep();    // move a vector field forward one time step
        void sStep();    // move a scalar field forward one time step

        inline int idx(int i, int j, int k){ return i + width * (j + height * k); }

    private:
        void setBoundary(float *u);
        void addForce(float *u, float *f, float dt, int flag);
        void addSource(float *u, float *source, float dt);
        void advect(float *u, float *u0, float **v, float dt);
        void diffuse(float *u, float *u0, float k, float dt);
        void project(float *u, float *u0, float dt);
        void dissipate(float *u, float *u0, float a, float dt);
}

#endif
