program karp
  ! karp.soln.f
  ! This simple program approximates pi by computing pi = integral
  ! from 0 to 1 of 4/(1+x*x)dx which is approximated by sum from
  ! k=1 to N of 4 / (1+((k-.5)/N)**2).  The only input data required is N.
  !
  ! 10/11/95 RLF  MPI Parallel version 1 
  ! 3/7/97   RLF  Replace nprocs and mynum with size and rank
  ! 3/21/97  RLF  Change floats to real*8
  ! SHM 8/29/97  Change input to read from a file to accommodate 
  !              VW Companion
  ! SHM 8/29/97  Replaced goto with do while
  ! Nils Smeds Aug 14, 2000  Converted to F90
  !
  ! Uses only the 6 basic MPI calls
  !
  implicit none
  integer, parameter:: DOUBLE=kind(1.0d0), SINGLE=kind(1.0)
  include "mpif.h"
  integer :: n,i,mpierr,rank,size,tag
  real(DOUBLE) :: err,pi,sum,w,x
  integer, dimension(MPI_STATUS_SIZE) :: status
  intrinsic atan

  pi = 4.0 * atan(1.0)
  tag = 111
  open (unit = 20,file = "values")
  !
  ! All processes call the startup routine to get their rank 
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD,size,mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierr)
  !
  ! -------  Each new approximation to pi begins here. -------------------
  ! (Step 1) Get first value of N 
  call solicit(n,size,rank)
  !
  ! (Step 2): do the computation in N steps
  ! Parallel Version: there are "size" processes participating.  Each
  ! process should do 1/size of the calculation.  Since we want
  ! i = 1..n but rank = 0, 1, 2..., we start off with rank+1.
  do while (n > 0)
    w = 1.0 / n
    sum = 0.0
    do i = rank+1,n,size
      sum = sum + f((i-0.5)*w)
    end do
    sum = sum * w
    !
    ! (Step 3): print the results  
    ! (Parallel version: collect partial results and let master process print it)
    if (rank == 0) then
      print *, "host calculated x=", sum
      do i = 1,size-1
        call MPI_Recv(x,1,MPI_DOUBLE_PRECISION,i,tag,MPI_COMM_WORLD,status, &
             mpierr)
        print *, "host got x=", x
        sum = sum + x
      end do
      err = sum - pi
      print *, "sum, err =", sum, err
    else
      call MPI_Send(sum,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,mpierr)
    end if
    ! Get a new value of N
    call solicit(n,size,rank)
  end do
  !
  call MPI_Finalize(mpierr)
  close (unit = 20)
  !
contains
  !
  real(DOUBLE) function f(x)
    implicit none
    real(DOUBLE), intent(in) :: x

    f = 4.0 / (1.0+x*x)
  end function f

  subroutine solicit(n,nprocs,mynum)
    ! Get a value for N, the number of intervals in the approximation
    ! (Parallel versions: master process reads in N and then
    ! sends N to all the other processes)
    ! Note: A single broadcast operation could be used instead, but
    ! is not one of the 6 basics calls.
    implicit none
    include "mpif.h"
    integer, intent(inout) :: n
    integer, intent(in) :: mynum,nprocs
    integer :: i,mpierr,tag
    integer, dimension(MPI_STATUS_SIZE) :: status
    tag = 112
    if (mynum == 0) then
      read (20,*) n
      print *, "Number of approximation intervals = ", n
      do i = 1,nprocs-1
        call MPI_Send(n,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,mpierr)
      end do
    else
      call MPI_Recv(n,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,mpierr) 
    end if
  end subroutine solicit
end program karp
