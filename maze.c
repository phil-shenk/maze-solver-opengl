#include "maze.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
 * print out maze to terminal
 */
void printCell(int k1, int k2, cell c){
	//print value
	int v = 0;
	if(k1==0)
		v += c.north;
	if(k1==3)
		v += c.south;
	if(k2==0)
		v += c.west;
	if(k2==3)
		v += c.east;

	// print char
	char p = ' ';
	if(v==1)
		p = '#';
	if(v==2)
		p = '#';
	printf("%c",p);

}

/*
 * print a maze
 */
void printMaze(cell maze[8][8]){
	int arraySize = 8;
	//printf("cell");
	for(int i=0; i<arraySize; i++){
		// 4 print rows for each maze row
		for(int k1=0; k1<4; k1++){ 
			for(int j=0; j<arraySize; j++){
				for (int k2=0; k2<4; k2++){
					cell c = maze[i][j];
					printCell(k1,k2,c);
				}
			}
			printf("\n");

		}
	}
}
void printVisited(int visited[8][8]){
	int arraySize = 8;
	for(int i=0; i<arraySize; i++){
		for(int j=0; j<arraySize; j++){
			printf("%d",visited[i][j]);
		}
		printf("\n");
	}
}

/*
 * populate a maze
 */
void generateMaze(cell maze[8][8]){
	int arraySize = 8;


	//seed random numbers (granularity = 1 second)
	srand(time(0));

	static int visited[8][8];

	//populate maze with all walls & unvisited to start
	for(int i=0; i<arraySize; i++){
		for(int j=0; j<arraySize; j++){
			maze[i][j].east = 1;
			maze[i][j].west = 1;
			maze[i][j].north= 1;
			maze[i][j].south= 1;
			visited[i][j] = 1;
		}
	}

	openCellRecursive(4,4,maze,visited);
	
}
/*
 * recursively open the next cell given the current cell
 * cr = current row
 * cc = current col
 * maze[][] = the maze of cell elements
 * visited[][] = flags to mark visited cells (1 unvisited, 0 visited)
 */
void openCellRecursive(int cr, int cc, cell maze[8][8], int visited[8][8]){
	
	//RECURSIVE MAZE GENERATOR ALGORITHM
	// 1 given current cell as param:
	// 2 mark current cell as visited
	visited[cr][cc] -= 1; // += for now to ensure none are visited twice.. 

	//printf("\n\nRECURSIVE CALL\n");
	//printMaze(maze);
	//printf("\nVISITED:\n");
	//printVisited(visited);

	// 3 while the current cell has any unvisited neighbor cells:
	int neighbors[4] = {0,0,0,0};
	// count free neighbors and populate neighbors[] arr
	int numFree = getNeighbors(cr, cc, visited, neighbors);
	while(numFree>0){
		// 	1 choose one of the unvisited neighbors
		// 	2 remove the wall between current cell and chosen cell
		int r = rand() % numFree;
		// must make sure the # random numbers are mapped correctly onto NESW
		for(; neighbors[r]==0; r++){
			//shouldnt need any content in this loop
		}
		int nr=cr;
		int nc=cc;
		if(r==0){
			//next is north
			nr -= 1;
			maze[cr][cc].north = 0;
			maze[nr][nc].south = 0;
		}
		if(r==1){
			//next is east
			nc += 1;
			maze[cr][cc].east  = 0;
			maze[nr][nc].west  = 0;
		}
		if(r==2){
			//next is south
			nr += 1;
			maze[cr][cc].south = 0;
			maze[nr][nc].north = 0;
		}
		if(r==3){
			//next is west
			nc -= 1;
			maze[cr][cc].west  = 0;
			maze[nr][nc].east  = 0;
		}	

		// 	3 invoke routine recursively for chosen cell
		openCellRecursive(nr, nc, maze, visited);
		// make sure to re-check neighbors after recursing i think
		neighbors[0]=0;
	       	neighbors[1]=0;
		neighbors[2]=0;
		neighbors[3]=0;
		numFree = getNeighbors(cr, cc, visited, neighbors);
	}
	return;
}

/*
 * count free neighbors and populate neighbors[] with 1's for free [N,E,S,W]
 */
int getNeighbors(int cr, int cc, int visited[8][8], int neighbors[4]){
	
	int numFreeNeighbors = 0;
	//check north
	if(cr>0){
		numFreeNeighbors += visited[cr-1][cc];
		neighbors[0]     += visited[cr-1][cc];
	}
	//check east
	if(cc<7){
		numFreeNeighbors += visited[cr][cc+1];
		neighbors[1]     += visited[cr][cc+1];
	}
	//check south
	if(cr<7){
		numFreeNeighbors += visited[cr+1][cc];
		neighbors[2]     += visited[cr+1][cc];
	}
	//check west
	if(cc>0){
		numFreeNeighbors += visited[cr][cc-1];
		neighbors[3]     += visited[cr][cc-1];
	}
	return numFreeNeighbors;

}


