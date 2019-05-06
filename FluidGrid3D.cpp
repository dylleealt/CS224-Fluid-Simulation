#include "FluidGrid.h"

FluidGrid::FluidGrid()
{
}

FluidGrid::~FluidGrid()
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

FluidGrid::init()
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

FluidGrid::reset()
{
    for (int i = 0; i < m_totalCells; ++i){
         m_vx[i] = 0.f;
         m_vy[i] = 0.f;
         m_vz[i] = 0.f;
         m_p[i] = 0.f;
         m_d[i] = 0.f;
    }
}
