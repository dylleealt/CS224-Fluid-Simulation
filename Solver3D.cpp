#include "Solver3D.h"

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
    delete p;
    delete p0;
    delete d;
    delete d0;
}

void Solver3D::init()
{
    numRows = 100;
    numCols = 100;
    numLayers = 100;
    totalSize = numRows * numCols * numLayers;
    visc = 0.3f;
    kS = 0.5f;
    aS = 0.3f;
    dt = 0.5f;

    vx = new float [totalSize];
    vy = new float [totalSize];
    vz = new float [totalSize];
    vx0 = new float [totalSize];
    vy0 = new float [totalSize];
    vz0 = new float [totalSize];

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

void Solver3D::reset()
{
    for (int i = 0; i < totalSize; ++i){
         vx[i] = 0.f;
         vy[i] = 0.f;
         vz[i] = 0.f;
         p[i] = 0.f;
         d[i] = 0.f;
    }
}
