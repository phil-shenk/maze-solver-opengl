#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "initShader.h"
#include "../linalglib/linalglib.h"

#include "maze.h"

#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))

// 2D vector struct to store UV texture coordinates
typedef struct
{
    GLfloat x;
    GLfloat y;
} vec2;

// total number of vertices
int num_vertices;

// mode controls the state of the animation
int mode = 0;

// animation counter
int count = 0;


// parameters for simulated person walking through maze
int walkerx = 0;
int walkerz = 0;
const float walkerHeight = 1.0;
const float walkerGaze = 0.5; 
int finalStep = 1;
int walkspeed = 100;

// transformation matrices
GLuint ctm_location;
mat4 ctm = {
	{1,0,0,0},
	{0,1,0,0},
	{0,0,1,0},
	{0,0,0,1}
};

GLuint model_view_location;
mat4 model_view_matrix = {
	{1,0,0,0},
	{0,1,0,0},
	{0,0,1,0},
	{0,0,0,1}
};

GLuint projection_location;
mat4 projection_matrix = {
	{1,0,0,0},
	{0,1,0,0},
	{0,0,1,0},
	{0,0,0,1}
};

int window_size = 1000;

// voyage = the full journey to explore the maze the 1st time
int voyage[64*2]; // worst case = each cell in 8x8 grid is traversed twice
int voyageIndex = 0;

// solution = the optimal solution through the maze
int solution[64]; // optimal solution never backtracks, so worst case length = each grid traversed once
int solutionIndex = 0;

/*
 * populate the supplied 'cube' array with the vertices of a box of the desired style
 * style: 0 = pole, 1,2,3,4 = N,E,S,W wall
 */
void makeBox(vec4 cube[36], int style, int i, int j){
	vec4 v1 = {1,0,1,1};
	vec4 v2 = {1,0,-1,1};
	vec4 v3 = {-1,0,-1,1};
	vec4 v4 = {-1,0,1,1};
	vec4 v5 = {1,1,1,1};
	vec4 v6 = {1,1,-1,1};
	vec4 v7 = {-1,1,-1,1};
	vec4 v8 = {-1,1,1,1};

	vec4 ordered[36] = {
		v2,v1,v6,
		v6,v1,v5,
		v3,v2,v7,
		v7,v2,v6,
		v4,v3,v8,
		v8,v3,v7,
		v1,v4,v5,
		v5,v4,v8,
		v6,v5,v7,
		v7,v5,v8,
		v4,v1,v3,
		v3,v1,v2
	};
	
	// populate the supplied cube array with the desired box parameters
	for(int k=0; k<36; k++){
		// transform all 36 vertices to the position of the wall
		if(style==0){
			mat4 t = translation(j-0.5, 0, i-0.5);
			cube[k] = multiply_m4v4(t,
					multiply_m4v4(scaling(0.15,1.2,0.15), ordered[k]));
		}
		//N
		if(style==1){
			mat4 t = translation(j, 0, i-0.5);
			cube[k] = multiply_m4v4(t,
					multiply_m4v4(scaling(0.5,1.0,0.1), ordered[k]));
		}
		//E
		if(style==2){
			mat4 t = translation(j+0.5, 0, i);
			cube[k] = multiply_m4v4(t,
					multiply_m4v4(scaling(0.1,1.0,0.5), ordered[k]));
		}
		//S
		if(style==3){
			mat4 t = translation(j, 0, i+0.5);
			cube[k] = multiply_m4v4(t,
					multiply_m4v4(scaling(0.5,1.0,0.1), ordered[k]));
		}
		//W
		if(style==4){
			mat4 t = translation(j-0.5, 0, i);
			cube[k] = multiply_m4v4(t,
					multiply_m4v4(scaling(0.1,1.0,0.5), ordered[k]));
		}
		if(style==5){
			mat4 t = translation(j-2,0.0,i-2);
			cube[k] = multiply_m4v4(t,
					multiply_m4v4(scaling(0.5,0.1,0.5), ordered[k]));
		}
	}
}

/*
 * move box values into vertices, colors, and tex_coords to transfer to the graphics card
 */
void copyCube(int wallCounter, int texture, vec4 cube[36], vec4 vertices[], vec4 colors[], vec2 tex_coords[]){
	for(int k=0; k<36; k++){
		// copy the 36 values over to the 'vertices' array
		vertices[k+(36*wallCounter)] = cube[k];

		// and populate the color values into the 'colors' array
		float r1 = (float)rand()/(float)(RAND_MAX);
		float g1 = (float)rand()/(float)(RAND_MAX);
		float b1 = (float)rand()/(float)(RAND_MAX);
		vec4 c1 = {r1, g1, b1, 1};
		colors[k+(36*wallCounter)] = c1;

		// and populate the tex_coords
		vec2 tc[6] = {
			{0,0.5},
			{0.5,0.5},
			{0,0},
			{0,0},
			{0.5,0.5},
			{0.5,0}
		};

		float xoffset = 0.0;
		float yoffset = 0.0;
		// stone texture for pillars
		if(texture == 0){
			xoffset = 0.5;
			yoffset = 0.0;
		}
		// brick texture for walls
		if(texture == 1){
			xoffset = 0.0;
			yoffset = 0.0;
		}
		// grass texture for floor
		if(texture == 2){
			xoffset = 0.0;
			yoffset = 0.5;
		}
		for(int l=0; l<6; l++){
			tc[l].x += xoffset;
			tc[l].y += yoffset;
		}

		tex_coords[k+(36*wallCounter)] = tc[k%6];
	}
}

/*
 * recursive function to solve the maze
 */
int recSolveMaze(int solution[64], int solutionCounter[1], int voyage[64*2], int voyageCounter[1], int visited[8][8], cell maze[8][8], int r, int c, int rf, int cf){

	//given row r and col c

	// mark as visited
	visited[r][c] += 1;

	// ret 1 if reached target
	if(r==rf && c==cf){
		return 1;
	}

	// check north
	if(r>0){
		if(visited[r-1][c] == 0 && maze[r][c].north == 0){
			voyage[voyageCounter[0]] = 1;
			voyageCounter[0]++;
			if(recSolveMaze(solution, solutionCounter, voyage, voyageCounter, visited, maze, r-1, c, rf, cf) == 1){
				// this was the right move, so log this (tail recursion so this list is reversed, which is perfect actually)
				solution[solutionCounter[0]] = 3;
				solutionCounter[0]++;
				return 1;
			}
			voyage[voyageCounter[0]] = 3;
			voyageCounter[0]++;
		}
	}
	// check east
	if(c<7){
		if(visited[r][c+1] == 0 && maze[r][c].east== 0){
			voyage[voyageCounter[0]] = 2;
			voyageCounter[0]++;
			if(recSolveMaze(solution, solutionCounter, voyage, voyageCounter, visited, maze, r, c+1, rf, cf) == 1){
				solution[solutionCounter[0]] = 4;
				solutionCounter[0]++;
				return 1;
			}
			voyage[voyageCounter[0]] = 4;
			voyageCounter[0]++;
		}
	}
	// check south
	if(r<7){
		if(visited[r+1][c] == 0 && maze[r][c].south == 0){
			voyage[voyageCounter[0]] = 3;
			voyageCounter[0]++;
			if(recSolveMaze(solution, solutionCounter, voyage, voyageCounter, visited, maze, r+1, c, rf, cf) == 1){
				solution[solutionCounter[0]] = 1;
				solutionCounter[0]++;
				return 1;
			}
			voyage[voyageCounter[0]] = 1;
			voyageCounter[0]++;
		}
	}
	// check west
	if(c>0){
		if(visited[r][c-1] == 0 && maze[r][c].west == 0){
			voyage[voyageCounter[0]] = 4;
			voyageCounter[0]++;
			if(recSolveMaze(solution, solutionCounter, voyage, voyageCounter, visited, maze, r, c-1, rf, cf) == 1){
				solution[solutionCounter[0]] = 2;
				solutionCounter[0]++;
				return 1;
			}
			voyage[voyageCounter[0]] = 2;
			voyageCounter[0]++;
		}
	}

	// return 0 if none are empty
	return 0;
}

/*
 * recursively solve maze
 * solution[64] is a list of directions that make up the shortest path from ri to rc
 *  ^ big brain moment: i add the positions on the tail end of the recursive call, so
 * this array will be filled in reverse order, easy for recursion AND makes sense for walking BACK thru maze :)
 * voyage[64] is a list of directions that make up the path as it was explored from ri to rc
 */
void solveMaze(int solution[64], int voyage[64*2], cell maze[8][8], int ri, int ci, int rf, int cf){

	int visited[8][8];
	for(int i=0; i<8; i++){
		for(int j=0; j<8; j++){
			visited[i][j] = 0;
		}
	}

	int solutionCounter[] = {0};
	int voyageCounter[] = {0};

	int o = recSolveMaze(solution, solutionCounter, voyage, voyageCounter, visited, maze, ri, ci, rf, cf);
	if(o == 1){
		printf("\nsolved maze\n");
	}
}



/* 
 * build the 3D walls of the maze using the maze[8][8] array that maze.c generated
 */
void buildMaze(cell maze[8][8], vec4 vertices[], vec4 colors[], vec2 tex_coords[]){
	solveMaze(solution, voyage, maze, 0, 0, 7, 7);
	printf("\nvoyage:\n");
	for(int i=0; i<64*2; i++){
		printf("%d",voyage[i]);
	}
	printf("\noptimal final->initial solution:\n");
	for(int i=0; i<64; i++){
		printf("%d",solution[i]);
	}

	// make entrance and exit
	maze[0][0].west = 0;
	maze[7][7].east = 0;

	// populate maze walls
	int wallCounter = 0;
	for(int i=0; i<8; i++){
		for(int j=0; j<8; j++){
			// check North and West for each cell
			if(maze[i][j].north == 1){
				vec4 cube[36];
				makeBox(cube,1,i,j);
				copyCube(wallCounter,1,cube,vertices,colors,tex_coords);
				wallCounter += 1;
			}
			if(maze[i][j].west == 1){
				vec4 cube[36];
				makeBox(cube,4,i,j);
				copyCube(wallCounter,1,cube,vertices,colors,tex_coords);
				wallCounter += 1;
			}
			//check East for last cell in each row
			if(j == 7){
				if(maze[i][j].east == 1){
					vec4 cube[36];
					makeBox(cube,2,i,j);
					copyCube(wallCounter,1,cube,vertices,colors,tex_coords);
					wallCounter += 1;
				}
			}
			//check South for each cell in last row
			if(i == 7){
				if(maze[i][j].south == 1){
					vec4 cube[36];
					makeBox(cube,3,i,j);
					copyCube(wallCounter,1,cube,vertices,colors,tex_coords);
					wallCounter += 1;
				}
			}
		}
	}

	// populate pillars
	for(int i=0; i<9; i++){
		for(int j=0; j<9; j++){
			vec4 cube[36];
			makeBox(cube,0,i,j);
			copyCube(wallCounter,0,cube,vertices,colors,tex_coords);
			wallCounter += 1;
		}
	}

	// populate grass floor
	for(int i=0; i<12; i++){
		for(int j=0; j<12; j++){
			vec4 cube[36];
			makeBox(cube,5,i,j);
			copyCube(wallCounter,2,cube,vertices,colors,tex_coords);
			wallCounter += 1;
		}
	}

	return;
}

/*
 * count the walls in the maze that maze.c generated
 * need wallCount to know how many vertices to allocate
 */
int countWalls(cell maze[8][8]){
	int wallCount = 0;
	for(int i=0; i<8; i++){
		for(int j=0; j<8; j++){
			wallCount += maze[i][j].north;
			wallCount += maze[i][j].west;
			if(j==7)
				wallCount += maze[i][j].east;
			if(i==7)
				wallCount += maze[i][j].south;
		}
	}
	return wallCount;
}

void init(void)
{
	// INITIAL MODEL VIEW
	vec4 eyePoint = {3.5, 1, 3.5, 1};
	vec4 atPoint  = {3.5, 0, 3.5, 1};
	vec4 upVector = {0,   0,  -1, 0};

	model_view_matrix = look_at(eyePoint, atPoint, upVector);
	print_m4(model_view_matrix);
	float w = 8.0;
	float h = 8.0;
	float near = 0.105;
	float far = 15.201;

	//todo: investigate these parameters, i think far and near might be switched
 	projection_matrix = frustum(-w, w, -h, h, far, near);
	print_m4(projection_matrix);

	// generate a maze
	static cell maze[8][8];
	generateMaze(maze);
	printMaze(maze);

	//calculate how many vertices will be necessary
	int numWalls = countWalls(maze) - 2; // we'll remove 2 walls for entry/exit
	int numWallVertices = numWalls*36;

	int numPillars = (8+1)*(8+1);
	int numPillarVertices = numPillars*36;

	int numFloors = 12*12;
	int numFloorVertices = numFloors*36;

	printf("%dwallverts, %dpillarverts",numWallVertices, numPillarVertices);
	num_vertices = numWallVertices + numPillarVertices + numFloorVertices;

	vec4 vertices[num_vertices];
	vec4 colors[num_vertices];

	vec2 tex_coords[num_vertices];

	// build the maze walls and pillars using cubes
	// and set texture coords for the walls and pillars
	buildMaze(maze, vertices, colors, tex_coords);

	// LOAD TEXTURE
	int width  = 800;
	int height = 800;
	GLubyte my_texels[width][height][3];

	FILE *fp = fopen("texture.raw", "r");

	fread(my_texels, width * height * 3, 1, fp);
	fclose(fp);

	GLuint program = initShader("vshader.glsl", "fshader.glsl");
	glUseProgram(program);

	GLuint mytex[1];
	glGenTextures(1, mytex);
	glBindTexture(GL_TEXTURE_2D, mytex[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, my_texels);
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

	int param;
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &param);

	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	GLuint buffer;
	glGenBuffers(1, &buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors) + sizeof(tex_coords), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices), sizeof(colors), colors);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors), sizeof(tex_coords), tex_coords);

	GLuint vPosition = glGetAttribLocation(program, "vPosition");
	glEnableVertexAttribArray(vPosition);
	glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));

	GLuint vColor = glGetAttribLocation(program, "vColor");

	glVertexAttribPointer(vColor, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0 + sizeof(vertices));

	GLuint vTexCoord = glGetAttribLocation(program, "vTexCoord");
	glEnableVertexAttribArray(vTexCoord);
	glVertexAttribPointer(vTexCoord, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0 + (sizeof(vertices) + sizeof(colors)));

	GLuint texture_location = glGetUniformLocation(program, "texture");
	glUniform1i(texture_location, 0);
	printf("texture_location: %i\n", texture_location);

	model_view_location = glGetUniformLocation(program, "model_view_matrix");
	projection_location = glGetUniformLocation(program, "projection_matrix");
	printf("model_view_location: %i\n", model_view_location);
	printf("projection_location: %i\n", projection_location);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.7, 0.7, 1.0, 1.0);
	glDepthRange(1,0);
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUniformMatrix4fv(model_view_location, 1, GL_FALSE, (GLfloat *) &model_view_matrix);
    glUniformMatrix4fv(projection_location, 1, GL_FALSE, (GLfloat *) &projection_matrix);

    glDrawArrays(GL_TRIANGLES, 0, num_vertices);
    glutSwapBuffers();
}

void keyboard(unsigned char key, int mousex, int mousey)
{
	if(key == 'q')
	{
		glutLeaveMainLoop();
	}

        //glutPostRedisplay();
	if(key == 'n')
	{
		mode ++;
		printf("entering mode %d\n",mode);
		if(mode == 5){
			// these assignments are only necessary if you use 'n' to skip
			//walkerx = 7;
			//walkerz = 7;
		}
		if(mode == 7){
			// these assignments are only necessary if you use 'n' to skip
			//walkerx = 7;
			//walkerz = 7;
		}

	}
	if(key == '0'){
		mode = 0;
	}
	if(key == '1'){
		mode = 1;
	}
	if(key == '2'){
		mode = 2;
	}
	if(key == '3'){
		mode = 3;
	}
	if(key == '4'){
		mode = 4;
	}

	// check mode and reset counter as necessary if animations are restarting
	count = 0;
}

void mouse(int button, int state, int x, int y)
{
		if(button == 3)
		{
			ctm = multiply_m4(ctm, scaling(1.02,1.02,1.02));
			glutPostRedisplay();
		}
		if(button == 4)
		{
			ctm = multiply_m4(ctm, scaling(1/1.02,1/1.02,1/1.02));
			glutPostRedisplay();
		}

}

void reshape(int width, int height)
{
    glViewport(0, 0, window_size, window_size);
}

//mode 0 and 1 starting points
vec4 m0e = {3.5, 15.0, 3.5, 1.0};
vec4 m0a = {3.5, 0.0, 3.5, 1.0};
vec4 m0u = {0.0, 0.0,-1.0, 0.0};

// mode 2 starting points
const float circleRadius = 7.0;
vec4 m2e = {3.5, walkerHeight+5.0, circleRadius+3.5, 1.0};
vec4 m2a = {3.5, 0.0, 3.5, 1.0};
vec4 m2u = {0.0, 1.0, 0.0, 0.0};

// mode 3 starting points
vec4 m3e = {-1.0, walkerHeight, 0.0, 1.0};
vec4 m3a = { 0.0, walkerHeight-walkerGaze, 0.0, 1.0};
vec4 m3u = { 0.0, 1.0, 0.0, 0.0};

// simple linear interpolation between points pi and pf with parametrization 0<=t<=1
vec4 philerp(vec4 pi, vec4 pf, float t){
	vec4 vp = subtract_v4(pf, pi);
	vec4 vt = scale_v4(t, vp);
	vec4 pt = add_v4(pi, vt);

	// the point t*100% of the way between pi and pf
	return pt;
}

/*
 * rotational interpolation (get from pi to pf around the unit circle in the best direction)
 */
vec4 roterp(vec4 pi, vec4 pf, float t){
	//
	printf("\nroterp:\n");
	print_v4(pi);
	print_v4(pf);
	float cosTheta = dot_v4(normalize_v4(pi), normalize_v4(pf));
	printf("%f",cosTheta);
	float theta = acos(cosTheta);
	printf("%f",theta);
	float dtheta = t*theta;
	printf("%f",dtheta);

	vec4 result = multiply_m4v4(y_rotation(dtheta), normalize_v4(pi));
	return result;
}
/*
 * directional interpolation to rotate character
 */
vec4 direrp(int dir, int nextdir, vec4 thisat, float t){

	int dirdiff = nextdir - dir;
	int middir = nextdir;
	//int maxdir = dir;
	//if(dir<nextdir){
	//	maxdir = nextdir;
	//}
	// set middir to the avg if applicable, otherwise to nextdir
	if(dirdiff == 2 || dirdiff == -2){
		middir = (dir+nextdir)/2;
	}
	//now we have dir->middir->nextdir
	// each should have a direct corresponding vector that will be relative to thisat
	vec4 nesw[] = {
		{0,0,0,0},
		{0,0,-1,0},
		{1,0,0,0},
		{0,0,1,0},
		{-1,0,0,0}
	};
	vec4 vdir = nesw[dir];
	vec4 vmid = nesw[middir];
	vec4 vnext= nesw[nextdir];

	if(t<0.5){
		float tt = t*2;
		//lerp between dir and middir
		return add_v4(thisat, philerp(vdir, vmid, tt));
	}else{
		float tt = (t-0.5)*2;
		return add_v4(thisat, philerp(vmid,vnext, tt));
	}
}

void idle(void){
	// looking straight down mode
	if(mode == 0){
		model_view_matrix = look_at(m0e, m0a, m0u);
	}
	// flying from top to one corner (transition to circling mode)
	if(mode == 1){
		int steps = 300;
		float t = count/((float)steps);

		vec4 e = philerp(m0e, m2e, t);
		vec4 a = philerp(m0a, m2a, t);
		vec4 u = philerp(m0u, m2u, t);

		model_view_matrix = look_at(e, a, u);

		count++;
		if(count>steps){
			count = 0;
			mode++;
		}
	}

	// circling overhead mode
	if(mode == 2){
		// map count to range [0,2pi)
		int steps = 600;

		float theta = count*2*3.141592653/((float)steps);

		vec4 eyePoint = {circleRadius*sin(theta)+3.5,walkerHeight+5,circleRadius*cos(theta)+3.5,1};
		vec4 atPoint  = {3.5,walkerHeight-walkerGaze,3.5,1};
		vec4 upVector = {0,1,0,0};

		model_view_matrix = look_at(eyePoint, atPoint, upVector);

		count++;
		//flyCounter %= steps;
		if(count > steps){
			count = 0;
			mode++;
		}
	}

	// flying down from above
	if(mode == 3){
		int steps = 300;
		float t = count/((float)steps);

		vec4 e = philerp(m2e, m3e, t);
		vec4 a = philerp(m2a, m3a, t);
		vec4 u = philerp(m2u, m3u, t);

		model_view_matrix = look_at(e, a, u);

		count++;
		if(count>steps){
			count = 0;
			mode++;
		}
	}

	// voyaging thru maze
	if(mode == 4){
		// animation steps, not voyage steps...
		int steps = walkspeed;

		//step by step move through
		//
		// for each step:
		if ( voyage[voyageIndex+1] != 0){
			int dir = voyage[voyageIndex];
			//#TODO: make this nice and elegant using % and // (/?) to get dx and dy directly from dir
			int dx = 0;
			int dz = 0;
			//north
			if(dir == 1){
				dx = 0;
				dz = -1;
			}//east
			if(dir == 2){
				dx = 1;
				dz = 0;
			}//south
			if(dir == 3){
				dx = 0;
				dz = 1;
			}//west
			if(dir == 4){
				dx = -1;
				dz = 0;
			}
			vec4 thispos = {walkerx, walkerHeight, walkerz, 1};
			vec4 thisat = {walkerx+dx, walkerHeight-walkerGaze, walkerz+dz, 1};
			vec4 nextat  = {walkerx+(2*dx), walkerHeight-walkerGaze, walkerz+(2*dz), 1};
			vec4 nextpos = {walkerx+dx,  walkerHeight, walkerz+dz,  1};

			if(count<(steps/2)){
				float t = count/((float)(steps/2));


				// move with lerp
				vec4 e = philerp(thispos, nextpos, t);
				vec4 a = philerp(thisat,  nextat,  t);
				vec4 u = {0,1,0,0};

				model_view_matrix = look_at(e, a, u);

			} else {
			//2nd half of animation, THE TURN
				float t = (count-(steps/2))/((float)(steps/2));

				int nextdir = voyage[voyageIndex+1];
				//deprecated: use rotational interpolation (went with directional instead)
				/*
				//#TODO: make this nice and elegant using % and // (/?) to get dx and dy directly from dir
				// rotate to new dir with roterp
				vec4 turnedat = {thisat.x+ndx, thisat.y, thisat.z+ndz, 1};
				vec4 ri = subtract_v4(nextat, thisat);
				vec4 rf = subtract_v4(turnedat, thisat); //nextpos..
				vec4 ra = roterp(normalize_v4(ri),  normalize_v4(rf),  t);
				print_v4(ra);
				vec4 u = {0,1,0,0};
				//print_v4(e);
				//printf("%f",t);
				vec4 next = add_v4(thispos, ra);
				next.y -= 0.8;
				model_view_matrix = look_at(nextpos, next, u);
				*/

				// interpolate between the current and the desired direction to orient the character
				vec4 a_interp;
				if(dir != nextdir){
					a_interp = nextat;
					a_interp = direrp(dir, nextdir, thisat, t);
				}else{
				a_interp = nextat;
				}
				vec4 u = {0,1,0,0};
				model_view_matrix = look_at(nextpos, a_interp,  u);


			}

			count++;
			if(count>steps){
				count = 0;
				voyageIndex++;
				walkerx += dx;
				walkerz += dz;
			}

		}else{
			voyageIndex = 0;
			count = 0;
			mode++;
			mode++;
		}
	}

	if(mode == 5){
		// step out of far side of maze
		int steps = 100;
		float t = (count/((float)(steps)));

		vec4 ei = { walkerx, walkerHeight, walkerz, 1.0};
		vec4 ef = { walkerx, walkerHeight, walkerz, 1.0};
		vec4 ai = { walkerx+0.0, walkerHeight-walkerGaze, walkerz+0.0, 1.0};
		vec4 af = { walkerx-1.0, walkerHeight-walkerGaze, walkerz+0.0, 1.0};
		
		vec4 e = philerp(ei, ef, t);
		vec4 a = philerp(ai, af, t);
		vec4 u = {0,1,0,0};

		model_view_matrix = look_at(e, a, u);
		
		count++;
		if(count>steps){
			count = 0;
			mode++;
		}
	}

	// turn around 180degrees
	if(mode == 6){
		int steps=200;
		count++;
		float t = count/((float)steps);


		vec4 e = { 8.0, walkerHeight, 7.0, 1.0};
		// use interpolated a instead //vec4 a = { 9.0, walkerHeight-walkerGaze, 7.0, 1.0};
		vec4 u = { 0.0, 1.0, 0.0, 0.0};

		int dir = 2;
		int nextdir = 4;
		vec4 edown = e;
		edown.y -= walkerGaze;

		vec4 a_interp = direrp(dir, nextdir, edown, t);
		model_view_matrix = look_at(e, a_interp,  u);


		if(count>steps){
			count = 0;
			mode++;
		}
	}
	if(mode == 7){
		// step back into maze (done with 180 rotation)
		int steps = 50;
		float t = (count/((float)(steps)));

			int dir = solution[0];
			//#TODO: make this nice and elegant using % and // (/?) to get dx and dy directly from dir
			int dx = 0;
			int dz = 0;
			//north
			if(dir == 1){
				dx = 0;
				dz = -1;
			}//east
			if(dir == 2){
				dx = 1;
				dz = 0;
			}//south
			if(dir == 3){
				dx = 0;
				dz = 1;
			}//west
			if(dir == 4){
				dx = -1;
				dz = 0;
			}

		vec4 ei = { 8.0, walkerHeight, 7.0, 1.0};
		vec4 ef = { 7.0, walkerHeight, 7.0, 1.0};
		vec4 ai = { 7.0, walkerHeight-walkerGaze, 7.0, 1.0};
		vec4 af = { 7.0+dx, walkerHeight-walkerGaze, 7.0+dz, 1.0};
		
		vec4 e = philerp(ei, ef, t);
		vec4 a = philerp(ai, af, t);
		vec4 u = {0,1,0,0};

		walkerx = 7;
		walkerz = 7;
		model_view_matrix = look_at(e, a, u);
		
		count++;
		if(count>steps){
			count = 0;
			mode++;
		}
	}

	// reached end of maze
	if(mode == 8){
		// animation steps, not voyage steps...
		int steps = walkspeed;

		//step by step move through
		//
		// for each step:
		if ( solution[solutionIndex+1] != 0){
			int dir = solution[solutionIndex];
			//printf("dirfinal%d\n",dir);
			//#TODO: make this nice and elegant using % and // (/?) to get dx and dy directly from dir
			int dx = 0;
			int dz = 0;
			//north
			if(dir == 1){
				dx = 0;
				dz = -1;
			}//east
			if(dir == 2){
				dx = 1;
				dz = 0;
			}//south
			if(dir == 3){
				dx = 0;
				dz = 1;
			}//west
			if(dir == 4){
				dx = -1;
				dz = 0;
			}
			vec4 thispos = {walkerx, walkerHeight, walkerz, 1};
			vec4 thisat  = {walkerx+dx, walkerHeight-walkerGaze, walkerz+dz, 1};
			vec4 nextat  = {walkerx+(2*dx), walkerHeight-walkerGaze, walkerz+(2*dz), 1};
			vec4 nextpos = {walkerx+dx,  walkerHeight, walkerz+dz,  1};

			if(count<(steps/2)){
				float t = count/((float)(steps/2));


				// move with lerp
				vec4 e = philerp(thispos, nextpos, t);
				vec4 a = philerp(thisat,  nextat,  t);
				vec4 u = {0,1,0,0};

				model_view_matrix = look_at(e, a, u);

			} else {
			//2nd half of animation, THE TURN
				float t = (count-(steps/2))/((float)(steps/2));

				int nextdir = solution[solutionIndex+1];
				//#TODO: make this nice and elegant using % and // (/?) to get dx and dy directly from dir

				vec4 a_interp;
				if(dir != nextdir){
					a_interp = direrp(dir, nextdir, thisat, t);
				}else{
					a_interp = nextat;
				}
				vec4 u = {0,1,0,0};
				model_view_matrix = look_at(nextpos, a_interp,  u);


			}

			count++;
			if(count>steps){
				count = 0;
				solutionIndex++;
				//printf("solutionIndex %d\n",solutionIndex);
				walkerx += dx;
				walkerz += dz;
			}

		}else{
			//solutionIndex = 0;
			count = 0;
			mode++;
		}
	}
	
	if(mode == 9){
		// step out of beginning
		int steps = 100;
		float t = (count/((float)(steps)));

			int dir = solution[solutionIndex];
			//printf("dirfinal%d\n",dir);
			//#TODO: make this nice and elegant using % and // (/?) to get dx and dy directly from dir
			int dx = 0;
			int dz = 0;
			//north
			if(dir == 1){
				dx = 0;
				dz = -1;
			}//east
			if(dir == 2){
				dx = 1;
				dz = 0;
			}//south
			if(dir == 3){
				dx = 0;
				dz = 1;
			}//west
			if(dir == 4){
				dx = -1;
				dz = 0;
			}

		walkerx = 0-dx;
		walkerz = 0-dz;
		//printf("%dwx%dwz",walkerx,walkerz);
		vec4 ei = { walkerx, walkerHeight, walkerz, 1.0};
		vec4 ef = { 0.0, walkerHeight, 0.0 , 1.0};
		vec4 ai = { 0.0, walkerHeight-walkerGaze, 0.0, 1.0};
		vec4 af = { -1.0, walkerHeight-walkerGaze, 0.0, 1.0};
		
		vec4 e = philerp(ei, ef, t);
		vec4 a = philerp(ai, af, t);
		vec4 u = {0,1,0,0};

		model_view_matrix = look_at(e, a, u);
		
		count++;
		if(count>steps){
			count = 0;
			mode++;
		}
	}

	glutPostRedisplay();

}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_size, window_size);
    glutCreateWindow("maze time");
    glewInit();
    init();
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    glutMainLoop();

    return 0;
}
