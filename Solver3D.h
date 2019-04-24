#ifndef __SOLVER3D_H__
#define __SOLVER3D_H__

class Solver3D
{
    public:
        Solver3D();
	~Solver3d();
        void init();
	void reset();

	void vStep(); // move a vector field forward one time step 
	void sStep(); // move a scalar field forward one time step

    private:
	int width;
	int height;
	int depth;
	int totalSize; // total number of cells
	float visc;   // viscosity
	float kS;     // diffusion constant
	float aS;     // dissipation rate
	float dt;     // time step

	// vector fields
	float *vx;    // velocity components
	float *vy;
	float *vz;
	float *vx0;
	float *vy0;
	float *vz0;

	// we may not need these
	float *v[3];    // current velocity 
	float *v0[3];   // old velocity
	
	// scalar fields
	float *p;     // pressure
	float *p0;
	float *d;     // density
	float *d0;

	inline int idx(int i, int j, int k){ return i + width * (j + height * k); }

	void addSource(int flag);
	void diffuse(float *f, float *f0, float k);
	void advect(float *f, float *f0, float *v);
	void project(float *f, float *f0);
	void dissipate(float *f, float *f0);


#endif 
