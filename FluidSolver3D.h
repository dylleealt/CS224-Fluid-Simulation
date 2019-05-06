#ifndef __FLUIDSOLVER_H__
#define __FLUIDSOLVER_H__

class FluidSolver
{
    public:
        FluidSolver();
        ~FluidSolver();
        void init();
        void reset();

        float getWidth(){ return m_width; }
        float getHeight(){ return m_height; }
        float getDepth(){ return m_depth; }
        int getNumRows(){ return m_numRows; }
        int getNumCols(){ return m_numCols; }
        int getNumLayers(){ return m_numLayers; }
        int getTotalCells(){ return m_totalCells; }
        float getViscosity(){ return m_visc; }
        float getDiffusionRate(){ return m_kS; }
        float getDissipationRate(){ return m_aS; }
        float getTimeStep(){ return m_dt; }

        float **getVelocity(){ return m_v; }
        float *getPressure(){ return m_p; }
        float *getDensity(){ return m_d; }


        void vStep();    // move a vector field forward one time step
        void sStep();    // move a scalar field forward one time step

        inline int idx(int i, int j, int k){ return i + width * (j + height * k); }

    private:
        // origin is assumed to be at (0,0,0)
        float m_width;
        float m_height;
        float m_depth;
        int m_numCols;     // num cells along width
        int m_numRows;     // num cells along height
        int m_numLayers;   // num cells along depth
        int m_totalCells;  // total number of cells
        float m_visc;      // viscosity
        float m_kS;        // diffusion constant
        float m_aS;        // dissipation rate
        float m_dt;        // time step

        // vector fields
        float *m_vx;       // velocity components
        float *m_vy;
        float *m_vz;
        float *m_vx0;
        float *m_vy0;
        float *m_vz0;

        float *m_v[NDIM];     // current velocity
        float *m_v0[NDIM];    // old velocity

        // scalar fields
        float *m_p;        // pressure
        float *m_p0;
        float *m_d;        // density
        float *m_d0;

        void setBoundary(float *u);
        void addForce(float *u, float *f, float dt, int flag);
        void addSource(float *u, float *s, float dt);
        void advect(float *u, float *u0, float **v, float dt);
        void diffuse(float *u, float *u0, float k, float dt);
        void project(float *u, float *u0, float dt);
        void dissipate(float *u, float *u0, float rate, float dt);
}

#endif
