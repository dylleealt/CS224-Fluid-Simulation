#ifndef __RENDERER_H__
#define __RENDERER_H__

#include <memory>
#include <vector>
#include <QImage>

#include "glm.hpp"

#include "FluidSolver.h"
#include "LevelSetSolver.h"

class Renderer
{
    public:
        Renderer(float timestep, int numFrames, int numCells, float size, float radius, float visc, float diff, float rate, float vorticity, int numParticles);
        Renderer(const Renderer &renderer) {}
        ~Renderer();
        void reset();

        int getNumCols() { return m_numCols; }
        int getNumRows() { return m_numRows; }
        int getNumLayers() { return m_numLayers; }
//        int geNumCells() { return m_numCells; }
//        float getWidth() { return m_width; }
//        float getHeight() { return m_height; }
//        float getDepth() { return m_depth; }
//        float getHx() { return m_hx; }
//        float getHy() { return m_hy; }
//        float getHz() { return m_hz; }
        float getTimestep() { return m_timestep; }
        int getNumFrames() { return m_numFrames; }

        inline int idx(int i, int j, int k) { return i + m_numCols * (j + m_numRows * k); }

        void simulateAndRender(int width, int height);
        float rayMarch(glm::vec3 rayOrigin, glm::vec3 rayDir);

    private:
        float m_timestep;
        int m_numFrames;
        int m_numCells;

        std::unique_ptr<FluidSolver> m_fSolver;
        std::unique_ptr<LevelSetSolver> m_lsSolver;

        glm::mat3x3 relativeMatrix(glm::vec3 lookVector);
        float sdBox(glm::vec3 p, glm::vec3 b);

//        static float* initializeDistanceField(int numCells, float size, float radius);
//        static float* initializeDensityField(int numCells, float size, float radius);


        int m_numCols;     // num cells along each dimension
        int m_numRows;
        int m_numLayers;
//        int m_nx;          // num of cells not including the boundaries
//        int m_ny;
//        int m_nz;
//        int m_numCells;    // total number of cells

//        // origin is assumed to be at (0,0,0)
//        float m_width;     // simulation space in each dimension
//        float m_height;
//        float m_depth;
//        float m_hx;        // size of voxel in each dimension
//        float m_hy;
//        float m_hz;
//        float m_minX, m_maxX; // spatial boundaries
//        float m_minY, m_maxY;
//        float m_minZ, m_maxZ;

//        // constants
//        float m_r = 0.0;		// particle radius

//        // scalar fields
//        float* m_phi0; // old level set function values
//        float* m_phi;  // current level set function values

//        std::vector<glm::vec3> m_particles;

//        void addSource(float *u, float *s, float dt);
};

#endif
