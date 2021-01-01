#include "maze.h"

int main(int argc, char **argv)
{
	static cell cells[8][8];
	generateMaze(cells);
	printMaze(cells);

}
