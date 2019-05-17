#ifndef __LEVELSETSOLVER_H__
#define __LEVELSETSOLVER_H__


class LevelSetSolver
{
    public:
        LevelSetSolver(int numCols, int numRows, int numLayers, float width, float height, float depth, float dt, int numParticles);
        ~LevelSetSolver();
        void reset();

        int getNumCols() { return m_numCols; }
        int getNumRows() { return m_numRows; }
        int getNumLayers() { return m_numLayers; }
        int geNumCells() { return m_numCells; }
        float getWidth() { return m_width; }
        float getHeight() { return m_height; }
        float getDepth() { return m_depth; }
        float getHx() { return m_hx; }
        float getHy() { return m_hy; }
        float getHz() { return m_hz; }
        float getTimeStep(){ return m_dt; }

        float **getVelocity() { return m_v; }
        float *getPressure() { return m_p; }
        float *getDensity() { return m_d; }

        void vStep();
        void sStep();

        inline int idx(int i, int j, int k) { return i + m_numCols * (j + m_numRows * k); }
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
        float m_dt;				// time step
		float m_r = 0.0;		// particle radius

        // scalar fields
		std::vector<float> m_phi;  // current level set function values
		std::vector<float> m_phi0; // old level set function values
		
		std::vector<glm::vec3> m_particles;

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
