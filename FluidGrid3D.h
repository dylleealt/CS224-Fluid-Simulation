#ifndef __FLUIDGRID_H__
#define __FLUIDGRID_H__

inline int idx(int i, int j, int k){ return i + width * (j + height * k); }

class FluidGrid
{
    public:
        FluidGrid();
        ~FluidGrid();
        void init();
        void reset();

        float getWidth(){ return width; }
        float getHeight(){ return height; }
        float getDepth(){ return depth; }
        int getX(){ return X; }
        int getY(){ return Y; }
        int getZ(){ return Z; }
        int getTotalSize(){ return totalSize; }
        float getViscosity(){ return visc; }
        float getDiffusionRate(){ return kS; }
        float getDissipationRate(){ return aS; }
        float getTimeStep(){ return dt; }

        float **getVelocity(){ return v; }
        float *getPressure(){ return p; }
        float *getDensity(){ return d; }

    private:
        float width;
        float height;
        float depth;
        int X;
        int Y;
        int Z;
        int totalSize; // total number of cells
        float visc;   // viscosity
        float kS;     // diffusion constant
        float aS;     // dissipation rate
        float dt;     // time step

        // vector fields
        float *vx;    // velocity components
        float *vy;
        float *vz;
        float *vx0;
        float *vy0;
        float *vz0;

        float *v[3];    // current velocity
        float *v0[3];   // old velocity

        // scalar fields
        float *p;     // pressure
        float *p0;
        float *d;     // density
        float *d0;
}

#endif
