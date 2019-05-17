#include "Particles.h"

Particles::Particles()
{
	//TODO: initialize
}

Particles::~Particles()
{
//    delete positions;
}

//void Particles::initialize_particles() {
//	int particle_index = 0;
//	for (double i = 0.f; i < width; i++) {
//		for (double j = 0.f; j < height; j++) {
//			for (double k = 0.f; k < depth; k++) {
//				particles[particle_index] = glm::vec3(i, j, k);
//				particle_index++;
//			}
//		}
//	}
//}

//glm::vec3 Particles::interpolate_velocity(glm::vec3 pos, float **v) {
//	double x0 = floor(pos.x);
//	double x1 = ceil(pos.x);
//	double y0 = floor(pos.y);
//	double y1 = ceil(pos.y);
//	double z0 = floor(pos.z);
//	double z1 = ceil(pos.z);

//	double xd = (pos.x - x0) / (x1 - x0);
//	double yd = (pos.y - y0) / (y1 - y0);
//	double zd = (pos.z - z0) / (z1 - z0);

//	glm::vec3 c00 = (glm::vec3(v[0][x0], v[1][y0], v[2][z0]) * (1.0 - xd)) + (glm::vec3(v[0][x1], v[1][y0], v[2][z0]) * xd);
//	glm::vec3 c01 = (glm::vec3(v[0][x0], v[1][y0], v[2][z1]) * (1.0 - xd)) + (glm::vec3(v[0][x1], v[1][y0], v[2][z1]) * xd);
//	glm::vec3 c10 = (glm::vec3(v[0][x0], v[1][y1], v[2][z0]) * (1.0 - xd)) + (glm::vec3(v[0][x1], v[1][y1], v[2][z0]) * xd);
//	glm::vec3 c11 = (glm::vec3(v[0][x0], v[1][y1], v[2][z1]) * (1.0 - xd)) + (glm::vec3(v[0][x1], v[1][y1], v[2][z1]) * xd);

//	glm::vec3 c0 = (c00 * (1.0 - yd)) + (c10 * yd);
//	glm::vec3 c1 = (c01 * (1.0 - yd)) + (c11 * yd);

//	return (c0 * (1.0 - zd)) + (c1 * zd);
//}

//void Particles::update_positions(float **v) {
//	for (int i = 0; i < numParticles(); i++) {
//		glm::vec3 position = particles[i];
//		particles[i] = position + interpolate_velocity(position, v);
//	}
//}
