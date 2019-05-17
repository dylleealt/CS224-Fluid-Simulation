#include "LevelSetSolver.h"

#include <cmath>
#include <algorithm>
#include <string.h>
#include <iostream>

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
    m_phi(new float[numCols * numRows * numLayers]),
    m_phi0(new float[numCols * numRows * numLayers]),
	
	m_particles(numParticles)
{
	//TODO: make sure particles are initialized to lie on fluid boundary
}

LevelSetSolver::~LevelSetSolver()
{
    delete m_phi0;
    delete m_phi;
}

//inline float phi_p(glm::vec3 x_p, glm::vec3 x) {
//	return sqrt(SQUARE(x[0] - x_p[0]) + SQUARE(x[1] - x_p[1]) + SQUARE(x[2] - x_p[2])) - m_r;
//}

//glm::vec4 LevelSetSolver::closest_particle(glm::vec3 x) {
//	glm::vec3 closest = x;
//	float dist = inf;
//	for (glm::vec3 p : m_particles) {
//		dist_l = abs(glm::distance(x, p));
//		if (dist_l < dist) {
//			closest = p;
//			dist = dist_l;
//		}
//	}
//	return glm::vec4(closest, phi_p(closest, x));
//}

//void LevelSetSolver::initializePhi() {
//	//for each point in the phi-grid, take phi_p for the closest particle
//	for (int x = 0; x < m_numCols; x++) {
//		for (int y = 0; y < m_numRows; y++) {
//			for (int z = 0; z < m_numLayers; z++) {
//				glm::vec4 closest_p = closest_particle(glm::vec3(x,y,z));
//				glm::vec3 local_closest_particle = closest_p.xyz;
//				float local_phi_p = closest_p.w;
//                m_phi0[idx(x,y,z)] = local_phi_p;
//				//TODO: there are edge cases where this is flat-out wrong
//			}
//		}
//	}
//}


float LevelSetSolver::phi_offset(float* phi, int x, int y, int z, int offset, char axis) {
	if (axis == 'x') {
        return phi[idx(x+offset,y,z)];
	} else if (axis == 'y') {
        return phi[idx(x,y+offset,z)];
	} else {
        return phi[idx(x,y,z+offset)];
	}
}

//grid cell index of phi, partial derivative of phi WRT axis, which axis is to be used, fluid velocity at grid cell in the given axis
float LevelSetSolver::phi_axis(int x, int y, int z, float delta, char axis, float vel) {
	if (vel == 0) {
		return 0.0;
	} else {
        float v1;
        float v2;
        float v3;
        float v4;
        float v5;
//        if (delta == 0.0) {
//            delta = 1.0;
//        }
		if (vel > 0) {
            v1 = (phi_offset(m_phi0, x, y, z, -2, axis) - phi_offset(m_phi0, x, y, z, -3, axis)) / delta;
            v2 = (phi_offset(m_phi0, x, y, z, -1, axis) - phi_offset(m_phi0, x, y, z, -2, axis)) / delta;
            v3 = (phi_offset(m_phi0, x, y, z, 0, axis) - phi_offset(m_phi0, x, y, z, -1, axis)) / delta;
            v4 = (phi_offset(m_phi0, x, y, z, 1, axis) - phi_offset(m_phi0, x, y, z, 0, axis)) / delta;
            v5 = (phi_offset(m_phi0, x, y, z, 2, axis) - phi_offset(m_phi0, x, y, z, 1, axis)) / delta;
		} else {
            v1 = (phi_offset(m_phi0, x, y, z, 3, axis) - phi_offset(m_phi0, x, y, z, 2, axis)) / delta;
            v2 = (phi_offset(m_phi0, x, y, z, 2, axis) - phi_offset(m_phi0, x, y, z, 1, axis)) / delta;
            v3 = (phi_offset(m_phi0, x, y, z, 1, axis) - phi_offset(m_phi0, x, y, z, 0, axis)) / delta;
            v4 = (phi_offset(m_phi0, x, y, z, 0, axis) - phi_offset(m_phi0, x, y, z, -1, axis)) / delta;
            v5 = (phi_offset(m_phi0, x, y, z, -1, axis) - phi_offset(m_phi0, x, y, z, -2, axis)) / delta;
		}
        float epsilon = 0.00001;
        float S1 = ((13.0 / 12.0) * (SQUARE(v1 - (2.0 * v2) + v3)) + (0.25 * SQUARE(v1 - (4.0 * v2) + (3.0 * v3))));
        float S2 = ((13.0 / 12.0) * (SQUARE(v2 - (2.0 * v3) + v4)) + (0.25 * SQUARE(v2 - v4)));
        float S3 = ((13.0 / 12.0) * (SQUARE(v3 - (2.0 * v4) + v5)) + (0.25 * SQUARE((3.0 * v3) - (4.0 * v4) + v5)));
        float a1 = 0.1 * (1.0 / SQUARE(epsilon + S1));
        float a2 = 0.6 * (1.0 / SQUARE(epsilon + S2));
        float a3 = 0.3 * (1.0 / SQUARE(epsilon + S3));

        float a_sum = a1 + a2 + a3;

        float w1 = a1 / a_sum;
        float w2 = a2 / a_sum;
        float w3 = a3 / a_sum;
	
        return vel * ((w1 * ((v1 / 3.0) - ((7.0 * v2) / 6.0) + ((11.0 * v3) / 6.0)))
                    + (w2 * (((-1.0 * v2) / 6.0) + ((5.0 * v3) / 6.0) + (v4 / 3.0)))
                    + (w3 * ((v3 / 3.0) + ((5.0 * v4) / 6.0) - (v5 / 6.0))));
	}
}
float LevelSetSolver::phi_t(int x, int y, int z, glm::vec3 del, glm::vec3 vel) {
    return -1.0 * (phi_axis(x, y, z, del[0], 'x', vel[0]) + phi_axis(x, y, z, del[1], 'y', vel[1]) + phi_axis(x, y, z, del[2], 'z', vel[2]));
}

glm::vec3 LevelSetSolver::phi0_gradient(int x, int y, int z) {
    float grad_x = (m_phi0[idx(x + 1, y, z)] - m_phi0[idx(x - 1, y, z)]) * 0.5;
    float grad_y = (m_phi0[idx(x, y + 1, z)] - m_phi0[idx(x, y - 1, z)]) * 0.5;
    float grad_z = (m_phi0[idx(x, y, z + 1)] - m_phi0[idx(x, y, z - 1)]) * 0.5;
//    std::cerr<<"grad_x: "<<grad_x<<"; grad_y: "<<grad_y<<"; grad_z: "<<grad_z<<"\n";
    if (grad_x + grad_y + grad_z != 0.0) {
        return glm::normalize(glm::vec3(grad_x, grad_y, grad_z));
    }
    return glm::vec3(grad_x, grad_y, grad_z);
}

void LevelSetSolver::update(float **velocity) {
	//swap phi and phi0 so that after a step forward, phi is the most recent one
	//TODO: make sure types are super consistent
    float* temp = m_phi0;
    m_phi0 = m_phi;
	m_phi = temp;
	//for each point in the phi-grid, get phi_t
    for (int x = 3; x < m_nx - 2; x++) {
        for (int y = 3; y < m_ny - 2; y++) {
            for (int z = 3; z < m_nz - 2; z++) {
                float *vel0 = velocity[0];
                float *vel1 = velocity[1];
                float *vel2 = velocity[2];
//                std::cerr<<phi0_gradient(x,y,z).z<<" ";
//                std::cerr<<vel0[idx(x,y,z)]<<" "<<vel1[idx(x,y,z)]<<vel2[idx(x,y,z)]<<"\n";
                m_phi[idx(x,y,z)] = phi_t(x, y, z, phi0_gradient(x,y,z), glm::vec3(vel0[idx(x,y,z)],vel1[idx(x,y,z)],vel2[idx(x,y,z)]));
//                std::cerr<<m_phi[idx(x,y,z)]<<" ";
			}
            std::cerr<<"\n";
		}
	}
	
}


void LevelSetSolver::reset()
{
    // should work in most architectures
    memset(m_phi0, 0, sizeof(float) * m_numCells);
    memset(m_phi, 0, sizeof(float) * m_numCells);
}

//takes world space coords
float LevelSetSolver::interpolate(float *u, float x, float y, float z)
{
    // clamp to boundaries
    x = std::min(std::max(x, m_minX), m_maxX);
    y = std::min(std::max(y, m_minY), m_maxY);
    z = std::min(std::max(z, m_minZ), m_maxZ);
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
