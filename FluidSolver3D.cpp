#include "FluidSolver.h"

#include <cmath>

#define SWAP(f, f0) {auto tmp=f; f=f0; f0=tmp;}

FluidSolver::FluidSolver()
{
}

FluidSolver::~FluidSolver()
{
    delete vx;
    delete vy;
    delete vz;
    delete vx0;
    delete vy0;
    delete vz0;
    delete p;
    delete p0;
    delete d;
    delete d0;
}

FluidSolver::init()
{
    m_width = 100;
    m_height = 100;
    m_depth = 100;
    m_totalCells = m_width * m_height * m_depth;
    m_visc = 0.3f;
    m_kS = 0.5f;
    m_aS = 0.3f;
    m_dt = 1.f;

    vx = new float[m_totalCells];
    vy = new float[m_totalCells];
    vz = new float[m_totalCells];
    vx0 = new float[m_totalCells];
    vy0 = new float[m_totalCells];
    vz0 = new float[m_totalCells];

    m_v[0] = m_vx;
    m_v[1] = m_vy;
    m_v[2] = m_vz;
    m_v0[0] = m_vx0;
    m_v0[1] = m_vy0;
    m_v0[2] = m_vz0;

    m_p = new float[m_totalCells];
    m_p0 = new float[m_totalCells];
    m_d = new float[m_totalCells];
    m_d0 = new float[m_totalCells];
}

FluidSolver::reset()
{
    for (int i = 0; i < m_totalCells; ++i){
         m_vx[i] = 0.f;
         m_vy[i] = 0.f;
         m_vz[i] = 0.f;
         m_p[i] = 0.f;
         m_d[i] = 0.f;
    }
}

void FluidSolver::vStep()
{
    addForce();
    SWAP();
    advect();
    SWAP();
    diffuse();
    SWAP();
    project();
}

void FluidSolver::sStep()
{
    addSource();
    SWAP();
    advect();
    SWAP();
    diffuse();
}

void FluidSolver::setBoundary(float *v)
{

}

void FluidSolver::addForce(float *v, float *f, float dt, int flag)
{
    switch (flag){
      	case 1: // only gravity
      	    for (int i = 0; i < totalCells; ++i){
      		      // may need to adjust gravitational accel later due to units
      	        v[i] += -9.8f * dt;
      	    }
      	    break;
      	case 2: // only external force
            for (int i = 0; i < totalCells; ++i){
                v[i] += f[i] * dt;
      	    }
      	    break;
      	case 3: // both gravity and external force
      	    for (int i = 0; i < totalCells; ++i){
                v[i] += (f[i] - 9.8f) * dt;
      	    }
      	    break;
      	default: // do nothing
    }
}

void FluidSolver::addSource(float *u, float *source, float dt)
{
    for (int i = 0; i < totalCells; ++i){
        u[i] += source[i] * dt;
    }
}

void FluidSolver::advect(float *u, float *u0, float **v, dt)
{
    int curIdx;
    float x, y, z;

    for (int i = 0; i < width; ++i){
        for (int j = 0; j < height; ++j){
            for (int k = 0; k < depth; ++k){
                x = i + 0.5f;
                y = j + 0.5f;
                z = k + 0.5f;
                curIdx = idx(i, j, k);
                curPos = idx(x, y, z); // *D
                // add interpolation here and segment time steps
                // trace particle
                prevX = x - v[0][curIdx] * dt;
                prevY = y - v[1][curIdx] * dt;
                prevZ = z - v[2][curIdx] * dt;
                // clamp to boundaries
                prevX = min(max(prevX, 0.f), width - 1);
                prevY = min(max(prevY, 0.f), height - 1);
                prevZ = min(max(prevZ, 0.f), depth - 1);
                // update field
                oldPos = idx(prevX, prevY, prevZ);
                // add interpolation here
                u[curIdx] = u0[oldPos];
            }
        }
    }
    // set boundary
}

void FluidSolver::linSolve(float *u, float *u0, float a, float c)
{
    int numIterations = 20;
    for (int t = 0; t < numIterations; ++t){
        for (int i = 1; i < width - 1; ++i){
            for (int j = 1; j < height - 1; ++j){
                for (int k = 1; k < depth - 1; ++k){
                    u[idx(i,j,k)] = (u0[idx(i,j,k)] + a * (
                    u[idx(i-1,j,k)] + u[idx(i+1,j,k)] +
                    u[idx(i,j-1,k)] + u[idx(i,j+1,k)] +
                    u[idx(i,j,k-1)] + u[idx(i,j,k+1)])) / c
                }
            }
        }
    setBoundary(u);
    }
}

void FluidSolver::diffuse(float *u, float *u0, float k, float dt)
{
    float a = k * dt * (width - 2) * (height - 2) * (depth - 2);
    float c = 1 + 4 * a;
    linSolve(u, u0, a, c);
}

void FluidSolver::project(float *u, float *u0, float dt)
{
    for (int i = 1; i < width - 1; ++i){
        for (int j = 1; j < height - 1; ++j){
            for (int k = 1; k < depth - 1; ++k){
                // missing scale factor for divergence
                div[idx(i,j,k)] = -0.5f * (ux[idx(i+1,j,k)] - ux[idx(i-1,j,k)]
                + uy[idx(i,j+1,k)] - uy[idx(i,j-1,k)]
                + uz[idx(i,j,k+1)] - uz[idx(i,j,k-1)];
                p[idx(i,j,k)] = 0.f;
            }
        }
    }
    setBoundary(div);
    setBoundary(p);

    // how do a and c terms work here
    linSolve(p, div, 1, 4);
    for (int i = 1; i < width - 1; ++i){
        for (int j = 1; j < height - 1; ++j){
            for (int k = 1; k < depth - 1; ++k){
                // account for scale factor here too
                ux[idx(i,j,k)] -= 0.5f * (p[idx(i+1,j,k)] - p[idx(i-1,j,k)]);
                uy[idx(i,j,k)] -= 0.5f * (p[idx(i,j+1,k)] - p[idx(i,j-1,k)]);
                uz[idx(i,j,k)] -= 0.5f * (p[idx(i,j,k+1)] - p[idx(i,j,k-1)]);
            }
        }
    }
    setBoundary(ux);
    setBoundary(uy);
    setBoundary(uz);
}
