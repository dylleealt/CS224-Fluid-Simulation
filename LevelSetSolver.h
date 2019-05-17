#ifndef __LEVELSETSOLVER_H__
#define __LEVELSETSOLVER_H__

#include <memory>
#include <vector>

#include "glm.hpp"

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
        float getNx() { return m_nx; }
        float getNy() { return m_ny; }
        float getNz() { return m_nz; }
        float getHx() { return m_hx; }
        float getHy() { return m_hy; }
        float getHz() { return m_hz; }
        float getTimeStep(){ return m_dt; }

        float interpolate(float *u, float x, float y, float z);
        void update(float **velocity);

        inline int idx(int i, int j, int k) { return i + m_numCols * (j + m_numRows * k); }

        // scalar fields
        float* m_phi0; // old level set function values
        float* m_phi;  // current level set function values

    private:
        float phi_offset(float* phi, int x, int y, int z, int offset, char axis);
        float phi_axis(int x, int y, int z, float delta, char axis, float vel);

        float phi_t(int x, int y, int z, glm::vec3 del, glm::vec3 vel);
        glm::vec3 phi0_gradient(int x, int y, int z);

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
		
		std::vector<glm::vec3> m_particles;
};

#endif
