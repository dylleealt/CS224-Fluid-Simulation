## CS2240 Final Project
Implementing a subset of [Efficient Simulation of Large Bodies of Water by Coupling Two and Three Dimensional Techniques](https://graphics.pixar.com/library/TwoDThreeDWaterSim/paper.pdf)

#### Other resources
  * [Stable Fluids](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf)  
  * [Visual Simulation of Smoke](http://graphics.ucsd.edu/~henrik/papers/smoke/smoke.pdf) for vorticity confinement

#### Main files/directories
  * FluidSolver.cpp/.h simulates the velocity (and density field) of the fluid. Includes Navier-Stokes implementation from Stable Fluids and vorticity confinement.
  * LevelSetSolver.cpp/.h tracks the surface of the fluid
  * Renderer.cpp/.h includes raymarching code to render the level set
  * particles provides code to step forward through the fluid solver and render a frame at that number of time steps. (See Bugs/Issues)
  * visualization is an attempt at a GPU implementation of this project
  * fluids.mp4 is a video of the particle-fluid simulation. A swirling external force is being applied in this video.

#### Bugs/Issues
  * The level set diverges from the correct distance field very quickly, breaks after a few frames
  * Frames can only be rendered once per simulation in the particles project code, cannot repeatedly re-draw the scene. So we can't step forward -> render -> step forward -> render -> ... Rather we have to step forward a certain number of steps to get the frame we want. Rendering particles is O(n^2) in the number of frames which is very, very slow. get_images(1-4).py is a horrible attempt to alleviate this.

#### Group Members
Dylan Tian  
Thomas Vandermillen
