#ifndef __FLUIDGRID_H__
#define __FLUIDGRID_H__

#define NDIM 3

inline int idx(int i, int j, int k){ return i + width * (j + height * k); }

class FluidGrid
{
    public:
        FluidGrid();
        ~FluidGrid();
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
}

#endif
