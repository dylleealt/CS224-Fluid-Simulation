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
      
        void addSource(int flag);
        void diffuse(float *f, float *f0, float k);
        void advect(float *f, float *f0, float *v);
        void project(float *f, float *f0);
        void dissipate(float *f, float *f0);
}

#endif
