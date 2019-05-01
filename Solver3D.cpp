#include "Solver3D.h"

#include <cmath>

#define SWAP(f, f0) {float *tmp=f; f=f0; f0=tmp;}

Solver3D::Solver3D()
{
}

Solver3D::~Solver3D()
{
    delete vx;
    delete vy;
    delete vz;
    delete vx0;
    delete vy0;
    delete vz0;
	
	delete vel;
	delete vel0;
	delete particles;
	
    delete p;
    delete p0;
    delete d;
    delete d0;
}

void Solver3D::init()
{
    width = 100;
    height = 100;
    depth = 100;
    totalSize = width * height * depth;
    visc = 0.3f;
    kS = 0.5f;
    aS = 0.3f;
    dt = 1.f;

    vx = new float [totalSize];
    vy = new float [totalSize];
    vz = new float [totalSize];
    vx0 = new float [totalSize];
    vy0 = new float [totalSize];
    vz0 = new float [totalSize];
	
	particles = new glm::vec3 [totalSize]; //TODO: new variable for numParticles

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;
    v0[0] = vx0;
    v0[1] = vy0;
    v0[2] = vz0;

    p = new float [totalSize];
    p0 = new float [totalSize];
    d = new float [totalSize];
    d0 = new float [totalSize];
}

void Solver3D::initialize_particles() {
	int particle_index = 0;
	for (double i = 0.f; i < width; i++) {
		for (double j = 0.f; j < height; j++) {
			for (double k = 0.f; k < depth; k++) {
				particles[particle_index] = glm::vec3(i, j, k);
				particle_index++;
			}
		}
	}
}

glm::vec3 Solver3D::interpolate_velocity(glm::vec3 pos) {
	double x0 = floor(pos.x);
	double x1 = ceil(pos.x);
	double y0 = floor(pos.y);
	double y1 = ceil(pos.y);
	double z0 = floor(pos.z);
	double z1 = ceil(pos.z);
	
	double xd = (pos.x - x0) / (x1 - x0);
	double yd = (pos.y - y0) / (y1 - y0);
	double zd = (pos.z - z0) / (z1 - z0);
	
	glm::vec3 c00 = (glm::vec3(vx[x0], vy[y0], vz[z0]) * (1.0 - xd)) + (glm::vec3(vx[x1], vy[y0], vz[z0]) * xd);
	glm::vec3 c01 = (glm::vec3(vx[x0], vy[y0], vz[z1]) * (1.0 - xd)) + (glm::vec3(vx[x1], vy[y0], vz[z1]) * xd);
	glm::vec3 c10 = (glm::vec3(vx[x0], vy[y1], vz[z0]) * (1.0 - xd)) + (glm::vec3(vx[x1], vy[y1], vz[z0]) * xd);
	glm::vec3 c11 = (glm::vec3(vx[x0], vy[y1], vz[z1]) * (1.0 - xd)) + (glm::vec3(vx[x1], vy[y1], vz[z1]) * xd);
	
	glm::vec3 c0 = (c00 * (1.0 - yd)) + (c10 * yd);
	glm::vec3 c1 = (c01 * (1.0 - yd)) + (c11 * yd);
	
	return (c0 * (1.0 - zd)) + (c1 * zd);
}

void Solver3D::reset()
{
    for (int i = 0; i < totalSize; ++i){
         vx[i] = 0.f;
         vy[i] = 0.f;
         vz[i] = 0.f;
		 
		 vel[i] = glm::vec3(0.f,0.f,0.f);
		 
         p[i] = 0.f;
         d[i] = 0.f;
    }
	initialize_particles();
}

void Solver3D::vStep()
{
    SWAP(vx, vx0);
    SWAP(vy, vy0);
    SWAP(vz, vz0);
    addForce(1);	
    SWAP(vx, vx0);
    SWAP(vy, vy0);
    SWAP(vz, vz0);
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

void Solver3D::sStep()
{
}

void Solver3D::addSource(int flag)
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
		float wL = px[
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
