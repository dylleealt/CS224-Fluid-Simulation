#include "LevelSetSolver.h"

#include <cmath>
#include <algorithm>

#define SWAP(u, u0) {auto tmp=u; u=u0; u0=tmp;}
#define MIX(a, x, y) ((1 - a) * (x) + (a) * y)
#define SQUARE(x) ((x) * (x))

LevelSetSolver::LevelSetSolver(int numCols, int numRows, int numLayers, float width, float height, float depth, float dt, int numParticles):
	m_numCols(numCols),
	m_numRows(numRows),
	m_numLayers(numLayers),
	m_nx(numCols - 1),
	m_ny(numRows - 1),
	m_nz(numLayers - 1),
	m_numCells(numCols * numRows * numLayers),
	
	m_width(width),
	m_height(height),
	m_depth(depth),
	m_hx(width/numCols),
	m_hy(height/numRows),
	m_hz(depth/numLayers),
	m_minX(0.5 * (width/numCols)),
	m_minY(0.5 * (height/numRows)),
	m_minZ(0.5 * (depth/numLayers)),
	m_maxX((numCols - 1.5) * (width/numCols)),
	m_maxY((numRows - 1.5) * (height/numRows)),
	m_maxZ((numLayers - 1.5) * (depth/numLayers)),
	
	m_dt(dt),
	
	m_phi(numCols * numRows * numLayers),
	m_phi0(numCols * numRows * numLayers),
	
	m_particles(numParticles)
{
	//TODO: make sure particles are initialized to lie on fluid boundary
}

LevelSetSolver::~LevelSetSolver()
{
}

inline float phi_p(glm::vec3 x_p, glm::vec3 x) {
	return sqrt(SQUARE(x[0] - x_p[0]) + SQUARE(x[1] - x_p[1]) + SQUARE(x[2] - x_p[2])) - m_r;
}

glm::vec4 LevelSetSolver::closest_particle(glm::vec3 x) {
	glm::vec3 closest = x;
	float dist = inf;
	for (glm::vec3 p : m_particles) {
		dist_l = abs(glm::distance(x, p));
		if (dist_l < dist) {
			closest = p;
			dist = dist_l;
		}
	}
	return glm::vec4(closest, phi_p(closest, x));
}

void LevelSetSolver::initializePhi() {
	//for each point in the phi-grid, take phi_p for the closest particle
	for (int x = 0; x < m_numCols; x++) {
		for (int y = 0; y < m_numRows; y++) {
			for (int z = 0; z < m_numLayers; z++) {
				glm::vec4 closest_p = closest_particle(glm::vec3(x,y,z));
				glm::vec3 local_closest_particle = closest_p.xyz;
				float local_phi_p = closest_p.w;
				m_phi0[x][y][z] = local_phi_p;
				//TODO: there are edge cases where this is flat-out wrong
			}
		}
	}
}

//TODO: phi_values --> m_phi or m_phi0
inline double phi_offset(int x, int y, int z, int offset, char axis) {
	if (axis == 'x') {
		return phi_values[x+offset][y][z];
	} else if (axis == 'y') {
		return phi_values[x][y+offset][z];
	} else {
		return phi_values[x][y][z+offset];
	}
}

//grid cell index of phi, partial derivative of phi WRT axis, which axis is to be used, fluid velocity at grid cell in the given axis
double LevelSetSolver::phi_axis(int x, int y, int z, double delta, char axis, double vel) {
	if (vel == 0) {
		return 0.0;
	} else {
		double v1;
		double v2;
		double v3;
		double v4;
		double v5;
		if (vel > 0) {
			v1 = (phi_offset(x, y, z, -2, axis) - phi_offset(x, y, z, -3, axis)) / delta;
			v2 = (phi_offset(x, y, z, -1, axis) - phi_offset(x, y, z, -2, axis)) / delta;
			v3 = (phi_offset(x, y, z, 0, axis) - phi_offset(x, y, z, -1, axis)) / delta;
			v4 = (phi_offset(x, y, z, 1, axis) - phi_offset(x, y, z, 0, axis)) / delta;
			v5 = (phi_offset(x, y, z, 2, axis) - phi_offset(x, y, z, 1, axis)) / delta;
		} else {
			v1 = (phi_offset(x, y, z, 3, axis) - phi_offset(x, y, z, 2, axis)) / delta;
			v2 = (phi_offset(x, y, z, 2, axis) - phi_offset(x, y, z, 1, axis)) / delta;
			v3 = (phi_offset(x, y, z, 1, axis) - phi_offset(x, y, z, 0, axis)) / delta;
			v4 = (phi_offset(x, y, z, 0, axis) - phi_offset(x, y, z, -1, axis)) / delta;
			v5 = (phi_offset(x, y, z, -1, axis) - phi_offset(x, y, z, -2, axis)) / delta;
		}
		double epsilon = 0.00001;
		double a1 = 0.1 * (1.0 / square(epsilon + S1(v1, v2, v3)));
		double a2 = 0.6 * (1.0 / square(epsilon + S2(v2, v3, v4)));
		double a3 = 0.3 * (1.0 / square(epsilon + S3(v3, v4, v5)));

		double a_sum = a1 + a2 + a3;

		double w1 = a1 / a_sum;
		double w2 = a2 / a_sum;
		double w3 = a3 / a_sum;
	
		return vel * (w1 * ((v1 / 3.0) - ((7.0 * v2) / 6.0) + ((11.0 * v3) / 6.0)))
				   + (w2 * (((-1.0 * v2) / 6.0) + ((5.0 * v3) / 6.0) + (v4 / 3.0)))
				   + (w3 * ((v3 / 3.0) + ((5.0 * v4) / 6.0) - (v5 / 6.0)));
	}
}
double LevelSetSolver::phi_t(int x, int y, int z, glm::vec3 del, glm::vec3 vel) {
	return -1.0 * (phi_axis(x, y, z, del[0], 'x', vel[0]) + phi_axis(x, y, z, del[1], 'y', vel[1]) + phi_axis(x, y, z, del[2], 'z', vel[2]));
}

glm::vec3 LevelSetSolver::phi0_gradient(int x, int y, int z) {
	double grad_x = (m_phi0[x + 1][y][z] - m_phi[x - 1][y][z]) * 0.5;
	double grad_y = (m_phi0[x][y + 1][z] - m_phi[x][y - 1][z]) * 0.5;
	double grad_z = (m_phi0[x][y][z + 1] - m_phi[x][y][z - 1]) * 0.5;
	return glm::normalize(glm::vec3(grad_x, grad_y, grad_z));
}

void LevelSetSolver::update() {
	//swap phi and phi0 so that after a step forward, phi is the most recent one
	//TODO: make sure types are super consistent
	double* temp = m_phi0;
	m_phi0 = m_phi;
	m_phi = temp;
	//for each point in the phi-grid, get phi_t
	for (int x = 0; x < m_numCols; x++) {
		for (int y = 0; y < m_numRows; y++) {
			for (int z = 0; z < m_numLayers; z++) {
				m_phi[x][y][z] = phi_t(x, y, z, phi0_gradient(x,y,z), VELOCITY[x][y][z]); //TODO: get velocity
			}
		}
	}
	
}


void LevelSetSolver::reset()
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

void LevelSetSolver::vStep()
{
}

void LevelSetSolver::sStep()
{
}

float LevelSetSolver::interpolate(float *u, float x, float y, float z)
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
              MIX(t, u[idx(i0, j0 + 1. k0)], u[idx(i0, j0 + 1, k0 + 1)])
            ),
            MIX(s,
              MIX(t, u[idx(i0 + 1, j0, k0)], u[idx(i0 + 1, j0, k0 + 1)]),
              MIX(t, u[idx(i0 + 1, j0 + 1, k0)], u[idx(i0 + 1, j0 + 1, k0 + 1)])
            )
          );
}

void LevelSetSolver::setBoundary(float *u, int b)
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
    u[idx(m_nx, m_ny, m_nz)] = (u[idx(m_nx - 1, m_ny, m_nz)] + u[idx(m_nx, m_ny - 1, m_nz0)] + u[idx(m_nx, m_ny, m_nz - 1)]) / 3.f;
}

void LevelSetSolver::addForce(float *u, float *f, float dt, int flag)
{
    switch (flag){
      	case 1: // only gravity
      	    for (int i = 0; i < m_numCells; ++i){
      		      // may need to adjust gravitational accel later due to units
      	        u[i] += -9.8f * dt;
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
      	default: // do nothing
    }
}

void LevelSetSolver::addSource(float *u, float *source, float dt)
{
    for (int i = 0; i < m_numCells; ++i){
        u[i] += source[i] * dt;
    }
}

void LevelSetSolver::advect(float *u, float *u0, float **v, float dt, int b)
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

void LevelSetSolver::linSolve(float *u, float *u0, float a, float c, int b)
{
    static int numIterations = 20;
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
    setBoundary(u, b);
    }
}

void LevelSetSolver::diffuse(float *u, float *u0, float k, float dt, int b)
{
    float a = k * dt * (m_nx - 1) * (m_ny - 1) * (m_nz - 1);
    float c = 1 + 4 * a;
    linSolve(u, u0, a, c, b);
}

void LevelSetSolver::project(float **v, float *p, float *div)
{
    float *vx = v[0], *vy = v[1], *vz = v[2];
    for (int i = 1; i < width - 1; ++i){
        for (int j = 1; j < height - 1; ++j){
            for (int k = 1; k < depth - 1; ++k){
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
    linSolve(p, div, 1, 4);
    for (int i = 1; i < width - 1; ++i){
        for (int j = 1; j < height - 1; ++j){
            for (int k = 1; k < depth - 1; ++k){
                // account for scale factor here too
                vx[idx(i,j,k)] -= 0.5f * (p[idx(i+1,j,k)] - p[idx(i-1,j,k)]);
                vy[idx(i,j,k)] -= 0.5f * (p[idx(i,j+1,k)] - p[idx(i,j-1,k)]);
                vz[idx(i,j,k)] -= 0.5f * (p[idx(i,j,k+1)] - p[idx(i,j,k-1)]);
            }
        }
    }
    setBoundary(ux, 1);
    setBoundary(uy, 2);
    setBoundary(uz, 3);
}
