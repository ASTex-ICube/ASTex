#version 330 core
layout (location = 0) in vec4 position;

void main()
{
    gl_Position = vec4(position.xy * 2.0 - 1.0, 1.0, 1.0);
    gl_PointSize = 6.*(position.w + 1.);// 12.0;
}
