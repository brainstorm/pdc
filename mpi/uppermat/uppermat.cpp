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
#include <string>
#include <cstdlib>

#include "mpi.h"

#define N 4

using namespace std;

int main(int argc, char **argv) {
  int i,j;

  int n, rank, size, type = 99;
  MPI::Status status;

  MPI::Init(argc, argv);
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();

  // Send buffer
  char* buf;
  // Distance buffer
  double *receive_buf;

  // This is what we actually have to send
  vector<string> sequences;
  // Array to receive the distances in
  double* distances;
  

  if (rank == 0) { // Master
  
    for (int i=0; i<8; i++)
      sequences.push_back("slkdfjassdflkajsdfaow");
    
    // Length of a null-terminated sequence
    int seq_len = sequences[0].size()+1;
    int nr_sequences = sequences.size();
    int start_mess[2];
    start_mess[0] = nr_sequences;
    start_mess[1] = seq_len;
    MPI::COMM_WORLD.Bcast(start_mess, 2, MPI_INT, 0);
    
    buf = (char*)malloc(seq_len*nr_sequences*sizeof(char));
    
    // Copy the sequences into the buffer
    // XXX: Error prone!
    for (int i=0; i<nr_sequences; i++)
      strcpy(&buf[i*seq_len], sequences[i].c_str());
    
    // Broadcast all the sequences
    MPI::COMM_WORLD.Bcast(buf, seq_len*nr_sequences, MPI_CHAR, 0);
    
    // Total number of distances
    int total_nr_dists = (nr_sequences*nr_sequences-nr_sequences)/2;
    // Number of distances to be calculated
    int local_nr_dists = total_nr_dists/size;

    // Calculate distances
    // ...
    double local_dist_vec[local_nr_dists];

    // Dummy numbers..
    for (int i=0; i<local_nr_dists; i++)
      local_dist_vec[i] = rank*local_nr_dists+i;

    // Collect distances from all processes
    receive_buf = (double*)malloc(total_nr_dists*sizeof(double));
		MPI::COMM_WORLD.Gather(local_dist_vec, local_nr_dists, MPI::DOUBLE, 
                            receive_buf, local_nr_dists, MPI::DOUBLE, 0);

    for (int i=0; i<total_nr_dists; i++)
      cout << receive_buf[i] << " ";
    cout << std::endl;
  } 
  else {  // Worker
    int start_mess[2];
    // Receive the number of sequences and the sequence length
    MPI::COMM_WORLD.Bcast(start_mess, 2, MPI_INT, 0);
    
    int nr_sequences = start_mess[0];
    // Remember: length of null-terminated sequences!
    int seq_len = start_mess[1];

    buf = (char*)malloc(seq_len*nr_sequences*sizeof(char));
    
    // Receive the actual sequences
    MPI::COMM_WORLD.Bcast(buf, nr_sequences*seq_len, MPI::CHAR, 0);

    // Copy data from buffer to vector
    for (int i=0; i<nr_sequences; i++) 
      sequences.push_back(&buf[i*seq_len]);
	  
    // Here it is assumed that the total number of distances is divisible by nr of processes
    // This is normally not the case? But can be changed later
    
    // Number of distances calculated by this process
    int total_nr_dists = (nr_sequences*nr_sequences-nr_sequences)/2;
    int local_nr_dists = total_nr_dists/size;
    
    // Do calculations..

    // Since this is process _rank_, it should calculate: 
    // distance dist_vec[rank*(nr_distances/size)] to dist_vec[(rank+1)*(nr_distances/size)-1]
    // What sequence pairs this corresponds to can be calculated,
    //  and will be calculated later today

    double local_dist_vec[local_nr_dists];

    // Dummy numbers..
    for (int i=0; i<local_nr_dists; i++)
      local_dist_vec[i] = rank*local_nr_dists+i;
    
    
    // The master process collects all the data  
		MPI::COMM_WORLD.Gather(local_dist_vec, local_nr_dists, MPI::DOUBLE, 
                            NULL, NULL, NULL, 0);
  }

  MPI::Finalize();
  return 0;
}
