#include "FluidSolver.h"

#include <cmath>
#include <algorithm>
#include <string.h>

#define SWAP(u, u0) {auto tmp=u; u=u0; u0=tmp;}
#define MIX(a, x, y) ((1 - a) * (x) + (a) * y)

FluidSolver::FluidSolver(int x, int y, int z, float width, float height, float depth, float visc, float diff, float rate, float dt):
    m_numCols(x),
    m_numRows(y),
    m_numLayers(z),
    m_nx(x - 1),
    m_ny(y - 1),
    m_nz(z - 1),
    m_numCells(x * y * z),
    m_width(width),
    m_height(height),
    m_depth(depth),
    m_hx(width/x),
    m_hy(height/y),
    m_hz(depth/z),
    m_minX(0.5*(width/x)),
    m_minY(0.5*(height/y)),
    m_minZ(0.5*(depth/z)),
    m_maxX((x -  1.5) * (width/x)),
    m_maxY((y -  1.5) * (height/y)),
    m_maxZ((z -  1.5) * (depth/z)),
    m_visc(visc),
    m_kS(diff),
    m_aS(rate),
    m_dt(dt),
    m_vx(new float[m_numCells]),
    m_vy(new float[m_numCells]),
    m_vz(new float[m_numCells]),
    m_vx0(new float[m_numCells]),
    m_vy0(new float[m_numCells]),
    m_vz0(new float[m_numCells]),
    m_cx(new float[m_numCells]),
    m_cy(new float[m_numCells]),
    m_cz(new float[m_numCells]),
    m_p(new float[m_numCells]),
    m_div(new float[m_numCells]),
    m_d(new float[m_numCells]),
    m_d0(new float[m_numCells])
{
    m_v[0] = m_vx;
    m_v[1] = m_vy;
    m_v[2] = m_vz;
    m_v0[0] = m_vx0;
    m_v0[1] = m_vy0;
    m_v0[2] = m_vz0;
}

FluidSolver::~FluidSolver()
{
    delete m_vx;
    delete m_vy;
    delete m_vz;
    delete m_vx0;
    delete m_vy0;
    delete m_vz0;
    delete m_p;
    delete m_d;
    delete m_d0;
}

void FluidSolver::reset()
{
    // should work in most architectures
    memset(m_vx, 0, sizeof(float) * m_numCells);
    memset(m_vy, 0, sizeof(float) * m_numCells);
    memset(m_vz, 0, sizeof(float) * m_numCells);
    memset(m_p, 0, sizeof(float) * m_numCells);
    memset(m_d, 0, sizeof(float) * m_numCells);

    // guaranteed to work but is slower
    /*
    for (int i = 0; i < m_numCells; ++i){
         m_vx[i] = 0.f;
         m_vy[i] = 0.f;
         m_vz[i] = 0.f;
         m_p[i] = 0.f;
         m_d[i] = 0.f;
    }
    */
}

void FluidSolver::update(float visc, float diff, float rate, float dt, int flag)
{
    // update velocity
    addForce(m_vz, nullptr, dt, flag);

//    diffuse(m_vx0, m_vx, visc, dt, 1);
//    diffuse(m_vy0, m_vy, visc, dt, 2);
//    diffuse(m_vz0, m_vz, visc, dt, 3);

//    project(m_v0, m_p, m_div);

    advect(m_vx, m_vx0, m_v0, dt, 1);
    advect(m_vy, m_vy0, m_v0, dt, 2);
    advect(m_vz, m_vz0, m_v0, dt, 3);

    project(m_v, m_p, m_div);

    // update density
    addSource(m_d, m_d0, dt);
    diffuse(m_d0, m_d, diff, dt, 0);
    advect(m_d, m_d0, m_v, dt, 0);

    // currently not implementing dissipation
}

//takes world space coords
float FluidSolver::interpolate(float *u, float x, float y, float z)
{
    // get back to grid coordinates
    float i = x / m_hx;
    float j = y / m_hy;
    float k = z / m_hz;
    // indices
    int i0 = floor(i);
    int j0 = floor(j);
    int k0 = floor(k);
    // weights
    float r = i - i0;
    float s = j - j0;
    float t = k - k0;
    return MIX(r,
            MIX(s,
              MIX(t, u[idx(i0, j0, k0)], u[idx(i0, j0, k0 + 1)]),
              MIX(t, u[idx(i0, j0 + 1, k0)], u[idx(i0, j0 + 1, k0 + 1)])
            ),
            MIX(s,
              MIX(t, u[idx(i0 + 1, j0, k0)], u[idx(i0 + 1, j0, k0 + 1)]),
              MIX(t, u[idx(i0 + 1, j0 + 1, k0)], u[idx(i0 + 1, j0 + 1, k0 + 1)])
            )
          );
}

void FluidSolver::setBoundary(float *u, int b)
{
    // along x-axis
    for (int y = 1; y < m_ny; ++y){
        for (int z = 1; z < m_nz; ++z){
            u[idx(0, y, z)] = (b == 1) ? -u[idx(1, y, z)] : u[idx(1, y, z)];
            u[idx(m_nx, y, z)] = (b == 1) ? -u[idx(m_nx - 1, y, z)] : u[idx(m_nx - 1, y, z)];
        }
    }
    // along y-axis
    for (int x = 1; x < m_nx; ++x){
        for (int z = 1; z < m_nz; ++z){
            u[idx(x, 0, z)] = (b == 2) ? -u[idx(x, 1, z)] : u[idx(x, 1, z)];
            u[idx(x, m_ny, z)] = (b == 2) ? -u[idx(x, m_ny - 1, z)] : u[idx(x, m_ny - 1, z)];
        }
    }
    // along z-axis
    for (int x = 1; x < m_nx; ++x){
        for (int y = 1; y < m_ny; ++y){
            u[idx(x, y, 0)] = (b == 3) ? -u[idx(x, y, 1)] : u[idx(x, y, 1)];
            u[idx(x, y, m_nz)] = (b == 3) ? -u[idx(x, y, m_nz - 1)] : u[idx(x, y, m_nz)];
        }
    }
    // corners
    u[idx(0, 0, 0)] = (u[idx(1, 0, 0)] + u[idx(0, 1, 0)] + u[idx(0, 0, 1)]) / 3;
    u[idx(m_nx, 0, 0)] = (u[idx(m_nx - 1, 0, 0)] + u[idx(m_nx, 1, 0)] + u[idx(m_nx, 0, 1)]) / 3.f;
    u[idx(0, m_ny, 0)] = (u[idx(1, m_ny, 0)] + u[idx(0, m_ny - 1, 0)] + u[idx(0, m_ny, 1)]) / 3.f;
    u[idx(0, 0, m_nz)] = (u[idx(1, 0, m_nz)] + u[idx(0, 1, m_nz)] + u[idx(0, 0, m_nz - 1)]) / 3.f;
    u[idx(m_nx, m_ny, 0)] = (u[idx(m_nx - 1, m_ny, 0)] + u[idx(m_nx, m_ny - 1, 0)] + u[idx(m_nx, m_nx, 1)]) / 3.f;
    u[idx(m_nx, 0, m_nz)] = (u[idx(m_nx - 1, 0, m_nz)] + u[idx(m_nx, 1, m_nz)] + u[idx(m_nx, 0, m_nz - 1)]) / 3.f;
    u[idx(0, m_ny, m_nz)] = (u[idx(1, m_ny, m_nz)] + u[idx(0, m_ny - 1, m_nz)] + u[idx(0, m_ny, m_nz - 1)]) / 3.f;
    u[idx(m_nx, m_ny, m_nz)] = (u[idx(m_nx - 1, m_ny, m_nz)] + u[idx(m_nx, m_ny - 1, m_nz)] + u[idx(m_nx, m_ny, m_nz - 1)]) / 3.f;
}

void FluidSolver::addForce(float *u, float *f, float dt, int flag)
{
    float G = -2.f;
    switch (flag){
      	case 1: // only gravity
      	    for (int i = 0; i < m_numCells; ++i){
      		      // may need to adjust gravitational accel later due to units
                u[i] += G * dt;
      	    }
      	    break;
      	case 2: // only external force
            for (int i = 0; i < m_numCells; ++i){
                u[i] += f[i] * dt;
      	    }
      	    break;
      	case 3: // both gravity and external force
      	    for (int i = 0; i < m_numCells; ++i){
                u[i] += (f[i] - 9.8f) * dt;
      	    }
      	    break;
        case 4: // swirl
            for (int i = 1; i < m_nx; ++i){
                for (int j = 1; j < m_ny; ++j){
                    for (int k = 1; k < m_nz; ++k){
                        int cx = m_nx / 2, cy = m_ny / 2;
                        float relx = i - cx, rely = j - cy;
                        float radius = relx * rely; // not really, we'll fix this later
                        // add an orthogonal vector to get swirl
                        m_vx[idx(i, j, k)] += -rely * dt / radius;
                        m_vy[idx(i, j, k)] += -relx * dt / radius;
                        m_vz[idx(i, j, k)] += G * dt;
                    }
                }
            }
//       	default: // do nothing
    }
}

void FluidSolver::addSource(float *u, float *source, float dt)
{
    for (int i = 0; i < m_numCells; ++i){
        u[i] += source[i] * dt;
    }
}

void FluidSolver::advect(float *u, float *u0, float **v, float dt, int b)
{
    int curIdx;
    float x, y, z, prevX, prevY, prevZ;

    for (int i = 1; i < m_nx; ++i){
        for (int j = 1; j < m_ny; ++j){
            for (int k = 1; k < m_nz; ++k){
                curIdx = idx(i, j, k);
                // get velocity at this point
                x = (i + 0.5f) * m_hx;
                y = (j + 0.5f) * m_hy;
                z = (k + 0.5f) * m_hz;
                // add multiple time steps
                // trace particle
                prevX = x - v[0][curIdx] * dt;
                prevY = y - v[1][curIdx] * dt;
                prevZ = z - v[2][curIdx] * dt;
                // clamp to boundaries
                prevX = std::min(std::max(prevX, m_minX), m_maxX);
                prevY = std::min(std::max(prevY, m_minY), m_maxY);
                prevZ = std::min(std::max(prevZ, m_minZ), m_maxZ);
                // update field
                u[curIdx] = interpolate(u0, prevX, prevY, prevZ);
            }
        }
    }
    // set boundary
    setBoundary(u, b);
}

void FluidSolver::linSolve(float *u, float *u0, float a, float c, int b)
{
    static int numIterations = 20;
    for (int t = 0; t < numIterations; ++t){
        for (int i = 1; i < m_nx; ++i){
            for (int j = 1; j < m_ny; ++j){
                for (int k = 1; k < m_nz; ++k){
                    u[idx(i,j,k)] = (u0[idx(i,j,k)] + a * (
                    u[idx(i-1,j,k)] + u[idx(i+1,j,k)] +
                    u[idx(i,j-1,k)] + u[idx(i,j+1,k)] +
                    u[idx(i,j,k-1)] + u[idx(i,j,k+1)])) / c;
                }
            }
        }
    setBoundary(u, b);
    }
}

void FluidSolver::diffuse(float *u, float *u0, float k, float dt, int b)
{
    float a = k * dt * (m_nx - 1) * (m_ny - 1) * (m_nz - 1);
    float c = 1 + 4 * a;
    linSolve(u, u0, a, c, b);
}

void FluidSolver::project(float **v, float *p, float *div)
{
    float *vx = v[0], *vy = v[1], *vz = v[2];
        for (int i = 1; i < m_nx; ++i){
            for (int j = 1; j < m_ny; ++j){
                for (int k = 1; k < m_nz; ++k){
                // missing scale factor for divergence
                div[idx(i,j,k)] = -0.5f * (vx[idx(i+1,j,k)] - vx[idx(i-1,j,k)]
                    + vy[idx(i,j+1,k)] - vy[idx(i,j-1,k)]
                    + vz[idx(i,j,k+1)] - vz[idx(i,j,k-1)]);
                p[idx(i,j,k)] = 0.f;
            }
        }
    }
    setBoundary(div, 0);
    setBoundary(p, 0);

    // how do a and c terms work here
    linSolve(p, div, 1, 4, 0);
    for (int i = 1; i < m_nx; ++i){
        for (int j = 1; j < m_ny; ++j){
            for (int k = 1; k < m_nz; ++k){
                // account for scale factor here too
                vx[idx(i,j,k)] -= 0.5f * (p[idx(i+1,j,k)] - p[idx(i-1,j,k)]);
                vy[idx(i,j,k)] -= 0.5f * (p[idx(i,j+1,k)] - p[idx(i,j-1,k)]);
                vz[idx(i,j,k)] -= 0.5f * (p[idx(i,j,k+1)] - p[idx(i,j,k-1)]);
            }
        }
    }
    setBoundary(vx, 1);
    setBoundary(vy, 2);
    setBoundary(vz, 3);
}
