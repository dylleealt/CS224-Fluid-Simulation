#include "FluidSolver.h"

#include <cmath>

#define SWAP(f, f0) {auto tmp=f; f=f0; f0=tmp;}

FluidSolver::FluidSolver()
{
}

FluidSolver::~FluidSolver()
{
}

void FluidSolver::vStep()
{
    SWAP(v, v0);
    addForce(1);
    SWAP(v, v0);
    advect(vx, vx0)
    advect(vy, vy0)
    advect(vz, vz0)
    SWAP(vx, vx0);
    SWAP(vy, vy0);
    SWAP(vz, vz0);
    diffuse(vx, vx0, visc);
    diffuse(vy, vy0, visc);
    diffuse(vz, vz0, visc);
    SWAP(vx, vx0);
    SWAP(vy, vy0);
    SWAP(vz, vz0);
    project(vx, vx0);
    project(vy, vy0);
    project(vz, vz0);
}

void FluidSolver::sStep()
{
}

void FluidSolver::addForce(float *fint flag)
{
    if (flag){
        // add gravity
	for (int i = 0; i < totalSize; ++i){
	    // may need to adjust gravitational constant later due to units
            vz[i] = vz0[i] + -9.8f * dt;
	}
	switch (flag){
            case 1: // add swirl
		int relx, rely;
		int cx = width / 2, cy = height / 2;
	        for (int i = 0; i < width; ++i){
                    for (int j = 0; j < height; ++j){
                       relx = i - cx;
		       rely = j - cy;
		       radius = relx * rely; // not really, we'll fix this later
		       // add an orthogonal vector to get swirl
                       vx[i] = vx0[i] + (-rely * dt / radius);
		       vy[i] = vy0[i] + (-relx * dt / radius);
		    }
		}
		break;
	}
    }
    // else no external forces
}

void advect(float *f, float *f0)
{
    float prevX, prevY, prevZ;
    for (int i = 1; i < width - 1; ++i){
        for (int j = 1; j < height - 1; ++j){
            for (int k = 1; k < width - 1; ++k){
		curIdx = idx(i, j, k);
                prevX = px[curIdx] - vx[curIdx] * dt;
                prevY = py[curIdx] - vy[curIdx] * dt;
                prevZ = pz[curIdx] - vz[curIdx] * dt;
		// clamp these values to the boundary
		int i0 = static_cast<int>(prevX - 0.5f);
		int j0 = static_cast<int>(prevY - 0.5f);
		int k0 = static_cast<int>(prevZ - 0.5f);
		int i1 = i0 + 1;
		int j1 = j0 + 1;
		int k1 = k0 + 1;
	    }
	}
    }
}

void diffuse(float *f, float *f0, float a)
{
    for (int i = 0; i < totalSize; ++i){
        f[i] = 0.f;
    }
    for (int t = 0; t < 10; ++t){
        for (int i = 1; i < width - 1; ++i){
            for (int j = 1; j < height; ++j){
	        for (int k = 1; k < depth; ++k){
                    f[idx(i,j,k)] = (f0[idx(i,j,k)] + a * (f[idx(i+1,j,k)] + f[idx(i-1,j,k)] + f[idx(i,j+1,k)] + f[idx(i,j-1,k)] + f[idx(i,j,k+1)] + f[idx(i,j,k-1)]) / (4.f * a + 1)

		}
	    }
	}
    }
}
