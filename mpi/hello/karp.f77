       program karp
c
c This simple program approximates pi by computing pi = integral
c from 0 to 1 of 4/(1+x*x)dx which is approximated by sum from
c k=1 to N of 4 / (1+((k-.5)/N)**2).  The only input data required is N.
c
c NOTE: Comments that begin with "cspmd" are hints for part b of the
c       lab exercise, where you convert this into an MPI program.
c
cspmd  Each process could be given a chunk of the interval to do.
c
c RLF 3/21/97  Change floats to real*8
c SHM 8/29/97  Change input to read from a file to accommodate 
c              VW Companion
c SHM 8/29/97  Replaced goto with do while
c
       real*8 err, f, pi, sum, w, x
       integer i, N
       f(x) = 4.0/(1.0+x*x)
       pi = 4.0*atan(1.0)
       open(unit=20, file='values')

cspmd  call startup routine that returns the number of tasks and the
cspmd  taskid of the current instance.

c Now read in a new value for N.  When it is 0, then you should depart.

       read(20,*) N
       print*, 'Number of approximation intervals = ',N

       do 20 while (N .GT. 0)
       w = 1.0/N
       sum = 0.0
       do i = 1,N
	  sum = sum + f((i-0.5)*w)
       enddo
       sum = sum * w
       err = sum - pi
       print *, 'sum = ',sum,' err =', err

       read(20,*) N
       print*, 'Number of approximation intervals = ',N

20     continue

       close(unit=20)
       end
