// parallel.cpp : Defines the entry point for the console application.




// sequence.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <assert.h>  
#include <string.h>
#define TAG 0

using namespace std;

const unsigned int GENS = 9;
unsigned int STARTROW=1,
	MAXROW = 9,
	MAXCOL = 9,
	GLOBAL_MAXROW = MAXROW,
	GLOBAL_MAXCOL = MAXCOL;


struct cell {
	int row;
	int col;
};

struct message_cell {
	int row;
	int col;
	int state;
};

enum cell_state { dead = 0, alive = 1 };

typedef std::vector<std::vector<cell_state>> Grid;
typedef std::vector<std::vector<int>> GridCount;
typedef std::vector<cell> List;
typedef std::vector<message_cell> MessageList;
typedef std::pair<int, int> ipair;
typedef std::vector<ipair> work_map;
typedef std::vector<std::vector<int>> int_int;
typedef void(*list_function)(cell);

Grid map;                // array holding cells
GridCount numNeighbors;  // array holding neighbor counts
List newlive,            // cells that have just been vivified
newdie,             // cells that have just died
maylive,            // candidates to vivify in next generation
maydie;           // candidates to kill in next generation
MessageList sendNewToTop,  // list contain new cells from first row
			sendNewToBottom; // list contain new cells from last row
MPI_Datatype  MPI_CELL;


void  MakeMpiCell(MPI_Datatype*  cell)
{
	// Set up cell datatype for communication
	MPI_Datatype    entryType[1] = { MPI_INT };
	int             count[1] = { 3 };
	MPI_Aint        addr[1] = { 0 };

	MPI_Type_struct(1, count, addr, entryType, cell);
	MPI_Type_commit(cell);
} // MakeMpiCell
	
void UpdateNeighbors(cell cc){
	int r,              // loop index for row of neighbor loops
		c;              // column loop index
	cell neighbor;      // structure form of a neighbor
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	cout << "Rank " << rank << " Updating neighbours for " << cc.row << " " << cc.col << std::endl;
	for (r = cc.row - 1; r <= cc.row + 1; ++r)
	for (c = cc.col - 1; c <= cc.col + 1; ++c){
		if ((r != cc.row || c != cc.col) && cc.row >= STARTROW && cc.row <= MAXROW &&
			cc.col >= 1 && cc.col <= MAXCOL) {  //skip cell itself and edge cases
			if (map[cc.row][cc.col] == alive){
				numNeighbors[r][c]++;
			}
			else {
				numNeighbors[r][c]--;
			}

			switch (numNeighbors[r][c]) {
			case 1:
				if (map[r][c] == alive) {
					neighbor.row = r;
					neighbor.col = c;
					maydie.push_back(neighbor);
				}
				break;
			case 3:
				if (map[r][c] == dead) {
					neighbor.row = r;
					neighbor.col = c;
					maylive.push_back(neighbor);
				}
				break;
			case 4:
				if (map[r][c] == alive) {
					neighbor.row = r;
					neighbor.col = c;
					maydie.push_back(neighbor);
				}
				break;
			}  // switch
		}
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
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
		cell c = { row, col };
		if (map[row][col] == dead) {
			newlive->push_back(c);
		}
		if (c.row == STARTROW){
			message_cell cñ = { c.row, c.col, (int)alive };
			sendNewToTop.push_back(cñ);
			cout << "Rank " << rank << " add cell to sendNewToTop " << c.row << " " << c.col << std::endl;
		}
		if (c.row == MAXROW){
			message_cell cñ = { c.row, c.col, (int)alive };
			sendNewToBottom.push_back(cñ);
			cout << "Rank " << rank << " add cell to sendNewToBottom " << c.row << " " << c.col << std::endl;
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
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::cout << "Rank: " << rank << " Generations: " << generations << std::endl;
	char* str = MapToString(map);
	std::cout << str;
}

void Vivify(cell c){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (map[c.row][c.col] == dead &&
		numNeighbors[c.row][c.col] == 3) {
		if (c.row >= STARTROW && c.row <= MAXROW &&
			c.col >= 1 && c.col <= MAXCOL) {
			map[c.row][c.col] = alive;
			newlive.push_back(c);
			if (c.row == STARTROW){
				message_cell cñ = { c.row, c.col, (int)alive };
				sendNewToTop.push_back(cñ);
				cout << "Rank " << rank << " add cell to sendNewToTop " << c.row << " " << c.col << std::endl;
			}
			if (c.row == MAXROW){
				message_cell cñ = { c.row, c.col, (int)alive };
				sendNewToBottom.push_back(cñ);
				cout << "Rank " << rank << " add cell to sendNewToBottom " << c.row << " " << c.col << std::endl;
			}
		}
	}
}

void Kill(cell c){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (map[c.row][c.col] == alive &&
		(numNeighbors[c.row][c.col] < 2 ||
		numNeighbors[c.row][c.col] > 3)) {
		if (c.row >= 1 && c.row <= MAXROW &&
			c.col >= 1 && c.col <= MAXCOL) {
			map[c.row][c.col] = dead;
			newdie.push_back(c);
			if (c.row == STARTROW){
				message_cell cc = { c.row, c.col, (int)dead };
				sendNewToTop.push_back(cc);
				cout << "Rank " << rank << " add cell to sendNewToTop "<< cc.row << " " << cc.col << std::endl;
			}
			if (c.row == MAXROW){
				message_cell cc = { c.row, c.col, (int)dead };
				sendNewToBottom.push_back(cc);
				cout << "Rank " << rank << " add cell to sendNewToBottom " << cc.row << " " << cc.col << std::endl;
			}
		}
	}
}

void Error(std::string str){
	std::cout << "Fatal: " << str << std::endl;
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

int myRankPlus1(){
	int rank, num_processes;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank + 1 < num_processes) {
		return rank + 1;
	}
	else {
		return MPI_PROC_NULL;
	}
}

int myRankMinus1(){
	int rank, num_processes;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank - 1 >= 0) {
		return rank - 1;
	}
	else {
		return MPI_PROC_NULL;
	}
}

void  UpdateBorders(message_cell* border, int size, Grid* Map)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	for (int i = 0; i < size; i++){
		message_cell mc = border[i];
		cell c = { mc.row, mc.col };
		cell_state state = (cell_state)mc.state;
		std::cout << "Rank " << rank << " add " << mc.row << " " << mc.col << " " << mc.state << std::endl;
		if (c.row >= STARTROW && c.row <= MAXROW &&
			c.col >= 1 && c.col <= MAXCOL)
			(*Map)[c.row][c.col] = state;
		if (state == alive) {
			//std::cout << "Rank " << rank << " Adding to newlive." << std::endl;
			newlive.push_back(c);
		}
		else {
			//std::cout << "Rank " << rank << " Adding to newdie." << std::endl;
			newdie.push_back(c);
		}
	}
	//std::cout << "Rank " << rank << " Finish update borders;" << std::endl;
} // UpdateBorders

void  SendToNeighbor(MessageList* sendListHigher,
					 MessageList* sendListLower,
					 Grid*        Map)
{
	int rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status    status1, status2;
	MPI_Request   sendHandleHigher;
	MPI_Request   recvHandleHigher;
	MPI_Request   sendHandleLower;
	MPI_Request   recvHandleLower;

	/*MessageList   recvListLower;
	MessageList   recvListHigher;
	recvListLower.resize(MAXCOL);
	recvListHigher.resize(MAXCOL);*/
	message_cell *sendCellsHigher,
		*sendCellsLower,
		*receiveCellsHigher = new message_cell[MAXCOL],
		*receiveCellsLower = new message_cell[MAXCOL];
	if (sendListHigher->size() == 0){
		sendListHigher->resize(1);
		sendListHigher->at(0).row = -2;
		sendListHigher->at(0).col = -2;
	}
	sendCellsHigher = &(*sendListHigher)[0];

	if (sendListLower->size() == 0) {
		sendListLower->resize(1);
		sendListLower->at(0).row = -2;
		sendListLower->at(0).col = -2;
	}
	sendCellsLower = &(*sendListLower)[0];

	int rankPlus1 = myRankPlus1(),
		rankMinus1 = myRankMinus1(),
		received;
	std::cout << "Rank" << rank << " Sending higher " << (*sendCellsHigher).row << " " << (*sendCellsHigher).col << std::endl;
	MPI_Isend(sendCellsHigher, sendListHigher->size(), MPI_CELL, rankMinus1, TAG, MPI_COMM_WORLD, &sendHandleHigher);
	std::cout << "Rank" << rank << " Sending lower " << (*sendCellsLower).row << " " << (*sendCellsLower).col << std::endl;
	MPI_Isend(sendCellsLower, 1, MPI_CELL, rankPlus1, TAG, MPI_COMM_WORLD, &sendHandleLower);

	MPI_Irecv(receiveCellsHigher, MAXCOL, MPI_CELL, rankMinus1, TAG, MPI_COMM_WORLD, &recvHandleHigher);
	MPI_Irecv(receiveCellsLower, MAXCOL, MPI_CELL, rankPlus1, TAG, MPI_COMM_WORLD, &recvHandleLower);

	MPI_Wait(&recvHandleHigher, &status1);
	MPI_Get_count(&status1, MPI_CELL, &received);
	//std::cout << "Rank" << rank << " Received higher " << received << std::endl;
	if (status1.MPI_SOURCE != MPI_PROC_NULL){
		//std::cout << std::endl;
		std::cout << "Rank" << rank << " Received " << receiveCellsHigher[0].row << " " << receiveCellsHigher[0].col << std::endl;
		UpdateBorders(receiveCellsHigher, received, Map);
		;
	}

	MPI_Wait(&recvHandleLower, &status2);
	MPI_Get_count(&status2, MPI_CELL, &received);
	//std::cout << "Rank" << rank << " Received lower " << received << std::endl;
	if (status2.MPI_SOURCE != MPI_PROC_NULL){
		//std::cout << std::endl;
		std::cout << "Rank" << rank << " Received " << receiveCellsLower[0].row << " " << receiveCellsLower[0].col << std::endl;
		UpdateBorders(receiveCellsLower, received, Map);
		;
	}
	
	MPI_Wait(&sendHandleHigher, &status2);
	MPI_Wait(&sendHandleLower, &status2);
} // SendToNeighbors

void simulate(){
	unsigned int gencount = 0;
	while (gencount < GENS) {
		WriteMap(map, gencount);
		gencount++;
		TraverseList(&maylive, Vivify);
		TraverseList(&maydie, Kill);
		maylive.clear();
		maydie.clear();
		SendToNeighbor(&sendNewToTop, &sendNewToBottom, &map);
		TraverseList(&newlive, UpdateNeighbors);
		TraverseList(&newdie, UpdateNeighbors);
		newlive.clear();
		newdie.clear();
		sendNewToBottom.clear();
		sendNewToTop.clear();
	}
	WriteMap(map, gencount);
}

int _tmain(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	int rank, num_processes, start, end, size;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int_int input;
	MakeMpiCell(&MPI_CELL);
	int* current_input;
	MPI_Status status;
	MPI_Request request;
	if (rank == 0) {
		std::string file_p("D:\\data.txt");
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
		MPI_Recv(&start, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&end, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&size, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
		//char * str = MapToString(map);
		//std::cout << "End:" << std::endl;
		//std::cout << str;
		for (int i = 1; i < num_processes; i++) {
			int str_size;
			char *str;
			MPI_Recv(&str_size, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			str = new char[str_size+1];
			//std::cout << "Receiving " << str_size << std::endl;
			MPI_Irecv(str, str_size, MPI_CHAR, i, TAG, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			std::cout << str;
			//delete[] str;
		}
	}
	else {
		char *buffer = MapToString(map);
		int str_size = strlen(buffer) + 1;
		//std::cout << "Sending: " << str_size << std::endl;
		MPI_Send(&str_size, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
		MPI_Isend(buffer, str_size, MPI_CHAR, 0, TAG, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
	}
	//
	MPI_Finalize();
}
