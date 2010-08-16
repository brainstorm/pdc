       program karp
c karp.soln.f
c This simple program approximates pi by computing pi = integral
c from 0 to 1 of 4/(1+x*x)dx which is approximated by sum from
c k=1 to N of 4 / (1+((k-.5)/N)**2).  The only input data required is N.
c
c 10/11/95 RLF  MPI Parallel version 1 
C 3/7/97   RLF  Replace nprocs and mynum with size and rank
c 3/21/97  RLF  Change floats to real*8
c SHM 8/29/97  Change input to read from a file to accommodate 
c              VW Companion
c SHM 8/29/97  Replaced goto with do while
c
c Uses only the 6 basic MPI calls
c
       real*8 pi, err, f, sum, w, x
       include 'mpif.h'
       integer i, N, rank, size
       integer status(MPI_STATUS_SIZE), mpierr, tag
       character*20 values
       f(x) = 4.0/(1.0+x*x)
       pi = 4.0*atan(1.0)
       tag = 111
       open(unit=20, file='values')

c All processes call the startup routine to get their rank 
       call MPI_INIT(mpierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierr)

c -------  Each new approximation to pi begins here. -------------------
c (Step 1) Get first value of N 
       call solicit (N,size,rank,'values')

c (Step 2): do the computation in N steps
c Parallel Version: there are "size" processes participating.  Each
c process should do 1/size of the calculation.  Since we want
c i = 1..n but rank = 0, 1, 2..., we start off with rank+1.
 
       do 20 while (N .GT. 0)
       w = 1.0/N
       sum = 0.0
       do i = rank+1,N,size
	  sum = sum + f((i-0.5)*w)
       enddo
       sum = sum * w

c (Step 3): print the results  
c (Parallel version: collect partial results and let master process print it)
       if (rank.eq.0) then
	  print *,'host calculated x=',sum
 	  do i = 1,size-1
 	     call MPI_RECV(x,1,MPI_DOUBLE_PRECISION,i,tag,
     .          MPI_COMM_WORLD,status,mpierr)
             print *,'host got x=',x
 	     sum=sum+x
 	  enddo
          err = sum - pi
          print *, 'sum, err =', sum, err
       else
 	  call MPI_SEND(sum,1,MPI_DOUBLE_PRECISION,0,tag,
     .        MPI_COMM_WORLD,mpierr)
       endif
c Get a new value of N
       call solicit (N,size,rank,'values')
20     continue       

       call MPI_FINALIZE(mpierr)
       close(unit=20)
       end

       subroutine solicit (N,nprocs,mynum,values)
c Get a value for N, the number of intervals in the approximation
c (Parallel versions: master process reads in N and then
c sends N to all the other processes)
c Note: A single broadcast operation could be used instead, but
c is not one of the 6 basics calls.
       include 'mpif.h'
       integer status(MPI_STATUS_SIZE), tag
       tag = 112
       if (mynum .eq. 0) then
          read(20,*) N
          print*, 'Number of approximation intervals = ',N
          do i = 1, nprocs-1
             call MPI_SEND(N,1,MPI_INTEGER,i,tag,
     .            MPI_COMM_WORLD,mpierr)
          enddo
       else
          CALL MPI_RECV(N,1,MPI_INTEGER,0,tag,
     .       MPI_COMM_WORLD,status,mpierr)
       endif
       return
       end
