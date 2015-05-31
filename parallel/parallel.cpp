// parallel.cpp : Defines the entry point for the console application.




// sequence.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <mpi.h>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <assert.h>  
#include <string.h>
#define TAG 0
const unsigned int GENS = 300;
unsigned int STARTROW=1,
	MAXROW = 9,
	MAXCOL = 9,
	GLOBAL_MAXROW = MAXROW,
	GLOBAL_MAXCOL = MAXCOL;


struct cell {
	int row;
	int col;
};

enum cell_state { dead = 0, alive = 1 };

typedef std::vector<std::vector<cell_state>> Grid;
typedef std::vector<std::vector<int>> GridCount;
typedef std::vector<cell> List;
typedef std::pair<int, int> ipair;
typedef std::vector<ipair> work_map;
typedef std::vector<std::vector<int>> int_int;
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


int get_worker_id_for(int row){
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int min_work = MAXROW / size;
	int remainder = MAXROW % size;
	int i = 0;
	if (remainder > 0) {
		min_work += 1;
		remainder -= 1;
	}

	while (row > min_work){
		int min_work = MAXROW / size;
		if (remainder > 0) {
			min_work += 1;
			remainder -= 1;
		}
		row -= min_work;
		i++;
	}
	return i;
}

void TraverseList(List* list, list_function f){
	for (int i = 0; i < list->size(); i++) {
		f(list->at(i));
	}
}

void ReadInput(std::string file_p, int_int *splittedTasks){
	int row = 0, col = 0,
		size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	splittedTasks->resize(size);
	std::ifstream file(file_p);

	file >> row;
	file >> col;
	//std::cout << row << "#" << col << std::endl;
	while (row != 0 && col != 0) {
		if (row >= 1 && row <= MAXROW && col >= 1 && col <= MAXCOL) {
			int id = get_worker_id_for(row);
			//std::cout << "Row " << row << " goes to worker " << id << std::endl;
			(*splittedTasks)[id].push_back(row);
			(*splittedTasks)[id].push_back(col);
		}
		file >> row >> col; // get next coordinate pair
	} // while
	file.close();
}

void Initialize(Grid& map, GridCount& neighbours, List* newlive, List* newdie, List* maylive, List* maydie, std::vector<int> *array) {
	int row = 0, col = 0;
	map.resize(MAXROW + 2);
	neighbours.resize(MAXROW + 2);
	for (int row = STARTROW-1; row <= MAXROW + 1; row++) {
		map[row].resize(MAXCOL + 2);
		neighbours[row].resize(MAXCOL + 2);
		for (int col = 0; col <= MAXCOL + 1; col++) {
			map[row][col] = dead;
			neighbours[row][col] = 0;
		}
	}

	for (int i = 0; i < array->size(); i += 2) {
		row = (*array)[i];
		col = (*array)[i+1];
		cell newCell = { row, col };
		if (map[row][col] == dead) {
			newlive->push_back(newCell);
		}
		map[row][col] = alive;
	}

	TraverseList(newlive, UpdateNeighbors);

	newlive->clear();
}

char* MapToString(Grid& map){
	int size = MAXCOL*(MAXROW - STARTROW + 1) + 1;
	char* str = new char [size];
	for (int i = 0; i < size; i++) str[i] = '\0';
	for (int row = STARTROW; row <= MAXROW; row++) {
		for (int col = 1; col <= MAXCOL; col++) {
			char* ch = (map[row][col] == alive) ? "A" : "_";
			str = strcat(str, ch);
		}
		str = strcat(str, "\n");
	}
	str = strcat(str, "\0");
	return str;
}

void WriteMap(Grid& map, int generations){
	std::cout << "Generations: " << generations << std::endl;
	char* str = MapToString(map);
	std::cout << str;
	delete[] str;
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

void simulate(){
	unsigned int gencount = 0;
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
}

work_map split_work(int size){
	work_map schedule;
	int min_work = MAXROW / size;
	int remainder = MAXROW % size;
	int next = 1;
	schedule.resize(size);
	for (int i = 0; i < size; i++){
		int end = next + min_work-1;
		if (remainder > 0) {
			end += 1;
			remainder -= 1;
		}
		ipair arr(next, end);
		schedule[i] = arr;
		next = end + 1;
	}
	return schedule;
}

void run(){
	int rank, num_processes, start, end, size;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int_int input;
	MPI_Datatype mpi_cell;
	int* current_input;

	if (rank == 0) {
		std::string file_p("Y:\\Documents\\data.txt");
		ReadInput(file_p, &input);
		current_input = (input[rank].size() > 0) ? &input[rank][0] : nullptr;
		size = input[rank].size();
		work_map schedule = split_work(num_processes);
		start = schedule[0].first;
		end = schedule[0].second;
		for (int dest_id = 1; dest_id < num_processes; dest_id++){
			int slave_start = schedule[dest_id].first,
				slave_end = schedule[dest_id].second;
			int* slave_input = &input[dest_id][0];
			//std::cout << "Frist element sent: " << slave_input[0] << std::endl;
			int size = input[dest_id].size();
			MPI_Send(&slave_start, 1, MPI_INT, dest_id, TAG, MPI_COMM_WORLD);
			MPI_Send(&slave_end, 1, MPI_INT, dest_id, TAG, MPI_COMM_WORLD);
			MPI_Send(&size, 1, MPI_INT, dest_id, TAG, MPI_COMM_WORLD);
			MPI_Send(slave_input, size, MPI_INT, dest_id, TAG, MPI_COMM_WORLD);
		}

	}
	else {
		MPI_Recv(&start, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&end, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&size, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		current_input = new int[size];
		STARTROW = start;
		MAXROW = end;
		//std::cout << "Received should be " << size << std::endl;
		MPI_Status status;
		int received;
		MPI_Recv(current_input, size, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &received);
		assert(received == size);
		/*std::cout << "Received ";
		for (int i = 0; i < size; i++){
			std::cout << current_input[i] << " ";
		}
		std::cout << std::endl;*/
	}
	STARTROW = start;
	MAXROW = end;
	std::vector<int> v(current_input, &current_input[size]);
	Initialize(map, numNeighbors, &newlive, &newdie, &maylive, &maydie, &v);
	simulate();
	if (rank == 0) {
		char * str = MapToString(map);
		std::cout << str;
		delete[] str;
		for (int i = 1; i < num_processes; i++) {
			int str_size;
			char *str;
			MPI_Recv(&str_size, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			str = new char[str_size];
			std::cout << "Receiving " << str_size << std::endl;
			MPI_Recv(str, str_size, MPI_CHAR, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::cout << str;
		}
	}
	else {
		char *buffer = MapToString(map);
		int str_size = strlen(buffer)+1;
				
		std::cout << "Sending: " << str_size << std::endl;
		MPI_Send(&str_size, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
		MPI_Send(buffer, str_size, MPI_CHAR, 0, TAG, MPI_COMM_WORLD);
		delete[] buffer;
	}
}

int _tmain(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	//clock_t begin = clock();
	//MPI_Barrier(MPI_COMM_WORLD);
	run();
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	MPI_Finalize();
	return 0;
}
