#version 400 core

in vec2 uv;
in vec3 color;

out vec4 fragColor;

void main() {
    // TODO [Task 19] Only color fragment if length(uv - 0.5) < 0.5, else discard
    if (length(uv - 0.5) < 0.5) {
       fragColor = vec4(color.xyz, 1);
    } else {
       discard;
    }
}
