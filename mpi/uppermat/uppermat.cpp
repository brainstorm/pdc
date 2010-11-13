// http://people.revoledu.com/kardi/tutorial/VB/Tips/Symmetric-Matrix.html
// http://www.thinkingparallel.com/2007/02/08/a-smart-way-to-send-a-stdvector-with-mpi-and-why-it-fails/
// http://www.mcs.anl.gov/research/projects/mpi/usingmpi/examples/simplempi/main.htm
// http://www.sgi.com/tech/stl/Multimap.html

// mpiCC uppermat.cpp && ./a.out

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <utility>

#include "mpi.h"

#define N 4

using namespace std;


void show(multimap<int, int> hash)
{

  multimap<int,int>::iterator it;
  pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ret;

  cout << "Contents of multimap (unrolled):" << endl;
  for (it = hash.begin(); it != hash.end(); ++it) {
	  cout << (*it).first << ":" << (*it).second << endl;
  }

  cout << "mymm contains:" << endl;
  for (int row=0; row < N; row++)
  {
    cout << row << " =>";
    ret = hash.equal_range(row);
    for (it=ret.first; it!=ret.second; ++it)
      cout << " " << (*it).second;
    cout << endl;
  }
}

int main(int argc, char **argv) {
  int i,j;
  multimap<int, int> rowhash;
  multimap<int, int>::iterator it;
  pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ret;

  int mat[N][N] = {{0, 1, 2, 3},
                   {0, 0, 4, 5},
                   {0, 0, 0, 6},
                   {0, 0, 0, 0}};


  int n, rank, size, type = 99;
  MPI::Status status;

  MPI::Init(argc, argv);
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();


  // Builds a vector with the upper matrix values, row-wise.
  // Does not include the diagonal line.
  for (i=0; i<N-1; i++) {
	  for (j=i+1; j<N; j++) {
		  rowhash.insert(pair<int, int>(i, mat[i][j]));
	  }
  }

  show(rowhash);
/*
  if (rank == 0) { //Master
	// sends rows to workers (ranks) so that:
	// rank = row = rowhash(key) ----> rowhash(values)
	for (i = rank + 1; i <= n; i += size) {
		  for (it = rowhash.begin(); it != rowhash.end(); ++it) {
			  cout << "Sending rows to workers...: " << (*it).first << ":" << (*it).second << endl;
		      MPI::COMM_WORLD.Send((*it).second, 1, MPI::INT, 0, 0);
		  }
	}
	// receives results from workers
	cout << n << endl;
  } else { // XXX: Rank matches the key (row nr) from the multimap, so each worker processes one row
	  MPI::COMM_WORLD.Reduce(&myrow, &n, rowhash.size, MPI::INT, MPI::SUM, 0);
	  // Sends result of the row to the master
	  MPI::COMM_WORLD.Send(&n, 1, MPI::INT, 0, 0);
  }
*/

  MPI::Finalize();
  return 0;
}
