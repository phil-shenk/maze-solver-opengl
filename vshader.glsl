#version 120

attribute vec4 vPosition;
attribute vec2 vTexCoord;

varying vec2 texCoord;

uniform mat4 model_view_matrix;
uniform mat4 projection_matrix;

void main()
{
	gl_Position = projection_matrix * model_view_matrix * vPosition/vPosition.w;
	texCoord = vTexCoord;

}
