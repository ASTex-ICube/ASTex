#version 330 core

#define MaxGraine 6

uniform vec2 uRes;

layout (std140) uniform uPoints
{
    vec4[MaxGraine] uGraines;
};

void main() {
    vec2 uv = (gl_FragCoord.xy/uRes.xy);

//    vec2[6] graines = vec2[6](vec2(0.81, 0.65), vec2(0.15, 0.25), vec2(0.52, 0.91), vec2(0.65, 0.28), vec2(0.12, 0.72), vec2(0.45, 0.48));
    vec3[6] colors = vec3[6](vec3(0., 0., 1.), vec3(0., 1., 0.), vec3(1., 0., 0.), vec3(0., 1., 1.), vec3(1., 0., 1.), vec3(1., 1., 0.));

    float dist = 1000.;
    vec3 color = vec3(1., 1., 1.);

    for(int g=0; g<MaxGraine; g++){
        vec2 graine = uGraines[g].xy;
        float wi = uGraines[g].z;
        float new_dist = length(uv-graine)*length(uv-graine) - wi;

//        if(new_dist == length(uv)){
//            break;
//        }

        if (new_dist <= dist) {
            color = colors[g%6];
            dist = new_dist;
        }
    }
    gl_FragColor = vec4(color, 1.0);
}
