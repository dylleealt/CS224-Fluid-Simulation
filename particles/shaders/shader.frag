#version 400 core

in vec3 WS_position; // world-space position
in vec3 WS_normal;   // world-space normal

out vec3 fragColor;

void main(){
    vec3 WS_toLight = normalize(vec3(10.0) - WS_position);
    float v = 0.2 + 0.5 * max(0.0, dot(normalize(WS_normal), WS_toLight));
    fragColor = vec3(0.5, 0.7, 0.956);
}
