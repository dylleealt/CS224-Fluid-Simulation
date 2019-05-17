#include "FluidSolver.h"

#include <cmath>
#include <algorithm>
#include <string.h>
#include <iostream>

#define SWAP(u, u0) {auto tmp=u; u=u0; u0=tmp;}
#define MIX(a, x, y) ((1 - a) * (x) + (a) * y)

FluidSolver::FluidSolver(int x, int y, int z, float width, float height, float depth, float visc, float diff, float rate, float vorticity, float dt) :
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
    m_maxZ((z - 1.5) * (depth/z)),
    m_visc(visc),
    m_kS(diff),
    m_aS(rate),
    m_swirl(vorticity),
    m_dt(dt)
{
    m_vx = new float[m_numCells]();
    m_vy = new float[m_numCells]();
    m_vz = new float[m_numCells]();
    m_vx0 = new float[m_numCells]();
    m_vy0 = new float[m_numCells]();
    m_vz0 = new float[m_numCells]();
    m_cx = new float[m_numCells]();
    m_cy = new float[m_numCells]();
    m_cz = new float[m_numCells]();
    m_vort = new float[m_numCells]();

    m_v[0] = m_vx;
    m_v[1] = m_vy;
    m_v[2] = m_vz;
    m_v0[0] = m_vx0;
    m_v0[1] = m_vy0;
    m_v0[2] = m_vz0;

    m_p = new float[m_numCells]();
    m_div = new float[m_numCells]();
    m_d = new float[m_numCells]();
    m_d0 = new float[m_numCells]();
    m_buf = new float[m_numCells]();
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

    memset(m_vx0, 0, sizeof(float) * m_numCells);
    memset(m_vy0, 0, sizeof(float) * m_numCells);
    memset(m_vz0, 0, sizeof(float) * m_numCells);
    memset(m_d0, 0, sizeof(float) * m_numCells);

    // guaranteed to work but is slower
    /*
    for (int i = 0; i < m_numCells; ++i){
         m_vx[i] = 0.f;
         m_vy[i] = 0.f;
         m_vz[i] = 0.f;
         m_p[i] = 0.f;
         m_d[i] = 0.f;
         m_vx0[i] = 0.f;
         m_vy0[i] = 0.f;
         m_vz0[i] = 0.f;
         m_d0[i] = 0.f;
    }
    */
}

void FluidSolver::update(float visc, float diff, float rate, float vorticity, float dt, int flag)
{
    // update velocity
    addForce(m_vx, m_vx0, dt, 2);
    addForce(m_vy, m_vy0, dt, 2);
    addForce(m_vz, m_vz0, dt, 4);

//    for (int i = 0; i < m_numCells; ++i){
//        std::cout<<"v: "<<m_vx[i]<<" "<<m_vy[i]<<" "<<m_vz[i]<<std::endl;
//    }

    vorticityConfinement(m_v, vorticity, dt);

//    for (int i = 0; i < m_numCells; ++i){
//        std::cout<<"v: "<<m_vx[i]<<" "<<m_vy[i]<<" "<<m_vz[i]<<std::endl;
//    }

    diffuse(m_vx0, m_vx, visc, dt, 1);
    diffuse(m_vy0, m_vy, visc, dt, 2);
    diffuse(m_vz0, m_vz, visc, dt, 3);

//    for (int i = 0; i < m_numCells; ++i){
//        std::cout<<"v: "<<m_vx0[i]<<" "<<m_vy0[i]<<" "<<m_vz0[i]<<std::endl;
//    }

    project(m_v0, m_p, m_div);

//    for (int i = 0; i < m_numCells; ++i){
//        std::cout<<"v: "<<m_vx0[i]<<" "<<m_vy0[i]<<" "<<m_vz0[i]<<std::endl;
//    }

    advect(m_vx, m_vx0, m_v0, dt, 1);
    advect(m_vy, m_vy0, m_v0, dt, 2);
    advect(m_vz, m_vz0, m_v0, dt, 3);

//    for (int i = 0; i < m_numCells; ++i){
//        std::cout<<"v: "<<m_vx[i]<<" "<<m_vy[i]<<" "<<m_vz[i]<<std::endl;
//    }

    project(m_v, m_p, m_div);

//    for (int i = 0; i < m_numCells; ++i){
//        std::cout<<"v: "<<m_vx[i]<<" "<<m_vy[i]<<" "<<m_vz[i]<<std::endl;
//    }

    // update density
    addSource(m_d, m_d0, dt);
    diffuse(m_d0, m_d, diff, dt, 0);
    advect(m_d, m_d0, m_v, dt, 0);

    // currently not implementing dissipation
}

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
    // std::cout<<i0<<" "<<j0<<" "<<k0<<std::endl;
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
    float G = -9.8;
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
                u[i] += (f[i] + G) * dt;
            }
            break;
        case 4: // gravity and swirl
            for (int i = 1; i < m_nx; ++i){
                for (int j = 1; j < m_ny; ++j){
                    for (int k = 1; k < m_nz; ++k){
                        int cx = m_nx / 2, cy = m_ny / 2;
                        float relx = i - cx, rely = j - cy;
                        float radius = std::max(1.f, relx + rely); // not really, we'll fix this later
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

// jacobi iteration solver
void FluidSolver::linSolve(float *u, float *u0, float a, float c, int b)
{
    static int numIterations = 80;
    for (int t = 0; t < numIterations; ++t){
        for (int i = 1; i < m_nx; ++i){
            for (int j = 1; j < m_ny; ++j){
                for (int k = 1; k < m_nz; ++k){
                    m_buf[idx(i,j,k)] = (u0[idx(i,j,k)] + a * (
                        u[idx(i-1,j,k)] + u[idx(i+1,j,k)] +
                        u[idx(i,j-1,k)] + u[idx(i,j+1,k)] +
                        u[idx(i,j,k-1)] + u[idx(i,j,k+1)])) / c;
                }
            }
        }
        SWAP(m_buf, u);
        setBoundary(u, b);
    }
}

void FluidSolver::diffuse(float *u, float *u0, float k, float dt, int b)
{
    float a = k * dt / (m_hx * m_hy * m_hz);
    float c = 1 + 6 * a;
    linSolve(u, u0, a, c, b);
}

void FluidSolver::project(float **v, float *p, float *div)
{
    float *vx = v[0], *vy = v[1], *vz = v[2];
        for (int i = 1; i < m_nx; ++i){
            for (int j = 1; j < m_ny; ++j){
                for (int k = 1; k < m_nz; ++k){
                    div[idx(i,j,k)] = -1.f/3.f * (
                        (vx[idx(i+1,j,k)] - vx[idx(i-1,j,k)]) / m_hx +
                        (vy[idx(i,j+1,k)] - vy[idx(i,j-1,k)]) / m_hy +
                        (vz[idx(i,j,k+1)] - vz[idx(i,j,k-1)]) / m_hz);
                    p[idx(i,j,k)] = 0.f;
            }
        }
    }
    setBoundary(div, 0);
    setBoundary(p, 0);

    linSolve(p, div, 1, 6, 0);
    for (int i = 1; i < m_nx; ++i){
        for (int j = 1; j < m_ny; ++j){
            for (int k = 1; k < m_nz; ++k){
                vx[idx(i,j,k)] -= 0.5f * (p[idx(i+1,j,k)] - p[idx(i-1,j,k)]) / m_hx;
                vy[idx(i,j,k)] -= 0.5f * (p[idx(i,j+1,k)] - p[idx(i,j-1,k)]) / m_hy;
                vz[idx(i,j,k)] -= 0.5f * (p[idx(i,j,k+1)] - p[idx(i,j,k-1)]) / m_hz;
            }
        }
    }
    setBoundary(vx, 1);
    setBoundary(vy, 2);
    setBoundary(vz, 3);
}

void FluidSolver::curl(float *c, float *u, float *v, int b)
{
    int du_dx = 0, du_dy = 0, du_dz = 0;
    int dv_dx = 0, dv_dy = 0, dv_dz = 0;

    // axes for curl
    if (b == 1) { du_dy = 1; dv_dz = 1; }
    else if (b == 2) { du_dz = 1; dv_dx = 1; }
    else if (b == 3) { du_dx = 1; dv_dy = 1; }

    for (int i = 1; i < m_nx; ++i){
        for (int j = 1; j < m_ny; ++j){
            for (int k = 1; k < m_nz; ++k){
                c[idx(i, j, k)] = 0.5f * (
                        (v[idx(i + du_dx, j + du_dy, k + du_dz)] - v[idx(i - du_dx, j - du_dy, k + du_dz)]) -
                        (u[idx(i + dv_dx, j + dv_dy, k + dv_dz)] - u[idx(i - dv_dx, j - dv_dy, k - dv_dz)])
                        );
            }
        }
    }
}

void FluidSolver::vorticityConfinement(float **v, float vorticity, float dt)
{
    float *vx = v[0], *vy = v[1], *vz = v[2];
    // calculate curl
    curl(m_cx, vy, vz, 1);
    curl(m_cy, vz, vx, 2);
    curl(m_cz, vx, vy, 3);

    // calculate magnitude of vorticity vector
    for (int i = 1; i < m_nx; ++i){
        for (int j = 1; j < m_ny; ++j){
            for (int k = 1; k < m_nz; ++k){
                int index = idx(i, j, k);
                m_vort[index] = sqrt(m_cx[index] * m_cx[index] + m_cy[index] * m_cy[index] + m_cz[index] * m_cz[index]);
            }
        }
    }

    setBoundary(m_vort, 0);

    // calculate gradient of vorticity
    for (int i = 1; i < m_nx; ++i){
        for (int j = 1; j < m_ny; ++j){
            for (int k = 1; k < m_nz; ++k){
                float gradvortx = 0.5f * (m_vort[idx(i + 1, j, k)] - m_vort[idx(i - 1, j, k)]);
                float gradvorty = 0.5f * (m_vort[idx(i, j + 1, k)] - m_vort[idx(i, j - 1, k)]);
                float gradvortz = 0.5f * (m_vort[idx(i, j, k + 1)] - m_vort[idx(i, j, k - 1)]);
                float norm = sqrt(gradvortx * gradvortx + gradvorty * gradvorty + gradvortz * gradvortz) + 0.0001f;
                // normalize gradient
                gradvortx /= norm;
                gradvorty /= norm;
                gradvortz /= norm;
                // add force
                vx[idx(i, j, k)] += vorticity * (gradvorty * m_cz[idx(i, j, k)] - gradvortz * m_cy[idx(i, j, k)]) * dt;
                vy[idx(i, j, k)] += vorticity * (gradvortz * m_cx[idx(i, j, k)] - gradvortx * m_cz[idx(i, j, k)]) * dt;
                vz[idx(i, j, k)] += vorticity * (gradvortx * m_cy[idx(i, j, k)] - gradvorty * m_cx[idx(i, j, k)]) * dt;
            }
        }
    }

    setBoundary(vx, 1);
    setBoundary(vy, 2);
    setBoundary(vz, 3);
}
