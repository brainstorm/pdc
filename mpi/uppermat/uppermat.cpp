// http://people.revoledu.com/kardi/tutorial/VB/Tips/Symmetric-Matrix.html
// http://www.thinkingparallel.com/2007/02/08/a-smart-way-to-send-a-stdvector-with-mpi-and-why-it-fails/
// http://www.mcs.anl.gov/research/projects/mpi/usingmpi/examples/simplempi/main.htm

#include <iostream>
#include <cmath>
#include <vector>
#include "mpi.h"

#define N 4

using namespace std;

int main(int argc, char **argv) {
  int i,j;
  vector<int> ary;

  int mat[N][N] = {{0, 2, 3, 2},
                   {0, 0, 1, 2},
                   {0, 0, 0, 2},
                   {0, 0, 0, 0}};


  int n, rank, size, type = 99;
  MPI::Status status;

  MPI::Init(argc, argv);
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();

  // Builds a vector with the upper matrix values, row-wise.
  // Does not include the matrix diagonal on the values.
  for (i=0; i<N-1; i++) {
	  for (j=i+1; j<N; j++) {
		  ary.push_back(mat[i][j]);
	  }
  }

 
  if (rank == 0) { //Master
	ary.resize(ary.size());
	MPI::COMM_WORLD.Recv(&n, sizeof(int), MPI::INT, 1, 0);
  } else {
	for (i=1; i<ary.size(); i++)
		MPI::COMM_WORLD.Send(&ary[0], i, MPI::INT, 0, 0);
  }

  MPI::Finalize();
  return 0;
}

/*
void showVec(vector<int> vec)
{
  for (vector<int>::iterator iter = vec.begin(); iter != vec.end(); ++iter) {
	  cout << *iter << endl;
  }
}
*/
