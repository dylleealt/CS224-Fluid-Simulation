#ifndef __FLUIDSOLVER_H__
#define __FLUIDSOLVER_H__


class FluidSolver
{
    public:
        FluidSolver();
        ~FluidSolver();
        void init();
        void reset();

        int getNumCols(){ return m_numCols; }
        int getNumRows(){ return m_numRows; }
        int getNumLayers(){ return m_numLayers; }
        int geNumCells(){ return m_numCells; }
        float getWidth(){ return m_width; }
        float getHeight(){ return m_height; }
        float getDepth(){ return m_depth; }
        float getHx() { return m_hx; }
        float getHy() { return m_hy; }
        float getHz() { return m_hz; }
        float getViscosity(){ return m_visc; }
        float getDiffusionRate(){ return m_kS; }
        float getDissipationRate(){ return m_aS; }
        float getTimeStep(){ return m_dt; }

        float **getVelocity(){ return m_v; }
        float *getPressure(){ return m_p; }
        float *getDensity(){ return m_d; }

        void vStep();
        void sStep();

        inline int idx(int i, int j, int k){ return i + m_numCols * (j + m_numRows * k); }
        float interpolate(float *u, float x, float y, float z);

    private:
        int m_numCols;     // num cells along each dimension
        int m_numRows;
        int m_numLayers;
        int m_nx;          // num of cells not including the boundaries
        int m_ny;
        int m_nz;
        int m_numCells;    // total number of cells

        // origin is assumed to be at (0,0,0)
        float m_width;     // simulation space in each dimension
        float m_height;
        float m_depth;
        float m_hx;        // size of voxel in each dimension
        float m_hy;
        float m_hz;
        float m_minX, m_maxX; // spatial boundaries
        float m_minY, m_maxY;
        float m_minZ, m_maxZ;

        // constants
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
        float *m_cx;       // curl compoments
        float *m_cy;
        float *m_cz;

        float *m_v[NDIM];  // current velocity
        float *m_v0[NDIM]; // old velocity

        // scalar fields
        float *m_p;        // pressure
        float *m_p0;
        float *m_d;        // density
        float *m_d0;

        void setBoundary(float *u, int b);
        void addForce(float *u, float *f, float dt, int flag);
        void addSource(float *u, float *s, float dt);
        void advect(float *u, float *u0, float **v, float dt, int b);
        void linSolve(float *u, float *u0, float a, float c, int b);
        void diffuse(float *u, float *u0, float k, float dt, int b);
        void project(float **v, float *p, float *div);
        void dissipate(float *u, float *u0, float rate, float dt);
}

#endif
