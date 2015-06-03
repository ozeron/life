// sequence.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>

const unsigned int MAXROW = 1000;
const unsigned int MAXCOL = 1000;
const unsigned int GENS = 300;

struct cell {
	int row;
	int col;
};

enum cell_state { dead = 0, alive = 1 };

typedef std::vector<std::vector<cell_state>> Grid;
typedef std::vector<std::vector<int>> GridCount;
typedef std::vector<cell> List;
typedef void(*list_function)(cell);

Grid map;                // array holding cells
GridCount numNeighbors;  // array holding neighbor counts
List newlive,            // cells that have just been vivified
	newdie,             // cells that have just died
	maylive,            // candidates to vivify in next generation
	maydie;             // candidates to kill in next generation


void UpdateNeighbors(cell c){
	int R,              // loop index for row of neighbor loops
		C;              // column loop index
	cell neighbor;      // structure form of a neighbor

	for (R = c.row - 1; R <= c.row + 1; ++R)
	for (C = c.col - 1; C <= c.col + 1; ++C)
	if (R != c.row || C != c.col) {  //skip cell itself
		if (map[c.row][c.col] == alive){
			numNeighbors[R][C]++;
		}
		else {
			numNeighbors[R][C]--;
		}

		switch (numNeighbors[R][C]) {
		case 1:
			if (map[R][C] == alive) {
				neighbor.row = R;
				neighbor.col = C;
				maydie.push_back(neighbor);
			}
			break;
		case 3:
			if (map[R][C] == dead) {
				neighbor.row = R;
				neighbor.col = C;
				maylive.push_back(neighbor);
			}
			break;
		case 4:
			if (map[R][C] == alive) {
				neighbor.row = R;
				neighbor.col = C;
				maydie.push_back(neighbor);
			}
			break;
		}  // switch
	}
}

void TraverseList(List* list, list_function f){
	for (int i = 0; i < list->size(); i++) {
		f(list->at(i));
	}
}

void Initialize(Grid& map, GridCount& neighbours, List* newlive, List* newdie, List* maylive, List* maydie, std::string file_p) {
	std::ifstream file(file_p);
	int row=0, col=0;
	map.resize(MAXROW+2);
	neighbours.resize(MAXROW+2);
	for (int row = 0; row <= MAXROW+1; row++) {
		map[row].resize(MAXCOL + 2);
		neighbours[row].resize(MAXCOL + 2);
		for (int col = 0; col <= MAXCOL+1; col++) {
			map[row][col] = dead;
			neighbours[row][col] = 0;
		}
	}

	file >> row;
	file >> col;
	//std::cout << row << "#" << col << std::endl;
	while (row != 0 && col != 0) {
		if (row >= 1 && row <= MAXROW && col >= 1 && col <= MAXCOL) {
			cell newCell = {row, col};
			if (map[row][col] == dead) {
				newlive->push_back(newCell);
			}
			map[row][col] = alive;
		}
		file >> row >> col; // get next coordinate pair
	} // while 
	file.close();

	TraverseList(newlive, UpdateNeighbors);

	newlive->clear();
}

void WriteMap(Grid& map, int generations){
	std::cout << "Generations: " << generations << std::endl;
	for (int row = 1; row <= MAXROW; row++) {
		for (int col = 1; col <= MAXCOL; col++) {
			char ch = (map[row][col] == alive) ? 'A' : '_';
			std::cout << ch;
		}
		std::cout << std::endl;
	}
}

void Vivify(cell c){
	if (map[c.row][c.col] == dead &&
		numNeighbors[c.row][c.col] == 3) {
		if (c.row >= 1 && c.row <= MAXROW &&
			c.col >= 1 && c.col <= MAXCOL) {
			map[c.row][c.col] = alive;
			newlive.push_back(c);
		}
	}	
}

void Kill(cell c){
	if (map[c.row][c.col] == alive &&
		(numNeighbors[c.row][c.col] < 2 || 
		 numNeighbors[c.row][c.col] > 3)) {
		if (c.row >= 1 && c.row <= MAXROW &&
			c.col >= 1 && c.col <= MAXCOL) {
			map[c.row][c.col] = dead;
			newdie.push_back(c);
		}
	}
}

void Error(std::string str){
	std::cout << "Fatal: " << str << std::endl;
}
int _tmain(int argc, _TCHAR* argv[]) {
	clock_t begin = clock();
	unsigned int gencount = 0;
	std::string file_p("Y:\\Documents\\data_b.txt");
	Initialize(map, numNeighbors, &newlive, &newdie, &maylive, &maydie, file_p);
	//WriteMap(map, gencount);
	while (gencount < GENS) {
		gencount++;
		TraverseList(&maylive, Vivify);
		TraverseList(&maydie, Kill);
		maylive.clear();
		maydie.clear();
		TraverseList(&newlive, UpdateNeighbors);
		TraverseList(&newdie, UpdateNeighbors);
		newlive.clear();
		newdie.clear();
	}
	//WriteMap(map, gencount);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Simulated " << GENS << " generations in " << elapsed_secs << " sec." << std::endl;
	return 0;
}