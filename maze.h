#ifndef _MAZE_H_
#define _MAZE_H_

typedef struct
{
	int north;
	int east;
	int south;
	int west;
} cell;

void printCell(int k1, int k2, cell c);
void printMaze(cell maze[8][8]);
void generateMaze(cell maze[8][8]);
void openCellRecursive(int cr, int cc, cell maze[8][8], int visited[8][8]);
int getNeighbors(int cr, int cc, int visited[8][8], int neighbors[4]);

#endif

