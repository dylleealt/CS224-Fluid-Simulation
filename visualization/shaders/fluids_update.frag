#version 400 core

uniform float firstPass;
uniform sampler2D prevPos;
uniform sampler2D prevVel;
uniform int gridSize;

// output from quad.vert
in vec2 uv;

// TODO [Task 15] setup the output locations
layout(location = 0) out vec4 pos;
layout(location = 1) out vec4 vel;

const float PI = 3.14159;
const float dt = 0.0167; // 1 sec/60 fps
const float G = -0.1;
const int N = 25;

/*
    A particle has the following properties:
        - pos.xyz = clip space position
        - pos.w = for now we don't use this
        - vel.xyz = clip space velocity
        - vel.w = pressure
*/

vec2 idx(int i, int j, int k){ // x, y, z indices
    return vec2(i + N * (j + N * k) / gridSize, 0);
}

// A helpful procedural "random" number generator
float hash(float n) { return fract(sin(n)*753.5453123); }

vec3 calculateInitialVelocity(int index) {
    float theta = PI * hash(index * 872.0238);
    const float MAX_VEL = 0.5;
    float velMag = MAX_VEL * hash(index * 98723.345);
    return velMag * vec3(cos(theta), sin(theta), cos(theta));
}

vec4 initPosition(int index) {
    return vec4(0);
}

vec4 initVelocity(int index) {
    return vec4(calculateInitialVelocity(index), 0);
}

vec4 updatePosition(int index) {
    vec4 curVel = texture(prevVel, uv);
    vec4 newPos = texture(prevPos, uv) + curVel * dt;
    newPos.w = 0;

    return newPos;
}

vec4 interpolate(vec4 p){
    // p /= h; // height of cell
    vec4 p0 = floor(p);
    int x0 = int(p0.x);
    int y0 = int(p0.y);
    int z0 = int(p0.z);
    // weights
    float r = p.x - p0.x;
    float s = p.y - p0.y;
    float t = p.z - p0.z;
    // neighboring cells
    return mix(
                mix(
                    mix(texture(prevVel, idx(x0, y0, z0)), texture(prevVel, idx(x0, y0, z0 + 1)), t),
                    mix(texture(prevVel, idx(x0, y0 + 1, z0)), texture(prevVel, idx(x0, y0 + 1, z0 + 1)), t),
                    s
                ),
                mix(
                    mix(texture(prevVel, idx(x0 + 1, y0, z0)), texture(prevVel, idx(x0 + 1, y0, z0 + 1)), t),
                    mix(texture(prevVel, idx(x0 + 1, y0 + 1, z0)), texture(prevVel, idx(x0 + 1, y0 + 1, z0 + 1)), t),
                    s),
                r
           );
}

vec4 addForce(vec4 p, vec4 v, int index, int flag) {
    return v + vec4(0, G, 0, 0) * dt;
}

vec4 advect(vec4 p, vec4 v, int index, int b) {
    p -= v * dt;
    // clamp to boundaries
    p = clamp(p, vec4(-1), vec4(1));
    v = interpolate(p);

    return v;
}

vec4 updateVelocity(int index) {
    const float G = -0.1;
    vec4 curPos = texture(prevPos, uv);
    vec4 curVel = texture(prevVel, uv);

    curVel = addForce(curPos, curVel, index, 0);
    curVel = advect(curPos, curVel, index, 0);

    return curVel;
}

void main() {
    int index = int(uv.x * gridSize);
    if (firstPass > 0.5) {
        pos = initPosition(index);
        vel = initVelocity(index);
    } else {
        pos = updatePosition(index);
        vel = updateVelocity(index);

//        if (pos.w < vel.w) {
//            pos = initPosition(index);
//            vel = initVelocity(index);
//        }
    }
}
