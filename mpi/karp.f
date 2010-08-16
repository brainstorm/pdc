program karp
  !
  ! This simple program approximates pi by computing pi = integral
  ! from 0 to 1 of 4/(1+x*x)dx which is approximated by sum from
  ! k=1 to N of 4 / (1+((k-.5)/N)**2).  The only input data required is N.
  !
  ! NOTE: Comments that begin with "cspmd" are hints for part b of the
  !       lab exercise, where you convert this into an MPI program.
  !
  !spmd  Each process could be given a chunk of the interval to do.
  !
  ! RLF 3/21/97  Change floats to real*8
  ! SHM 8/29/97  Change input to read from a file to accommodate 
  !              VW Companion
  ! SHM 8/29/97  Replaced goto with do while
  ! Nils Smeds Aug 14, 2000  Converted to F90
  !
  implicit none
  integer, parameter:: DOUBLE=kind(1.0d0), SINGLE=kind(1.0)
  integer :: n,i
  real(DOUBLE) :: err,pi,sum,w,x
  intrinsic atan

  pi = 4.0 * atan(1.0)
  open (unit = 20,file = "values")
  !
  !spmd  call startup routine that returns the number of tasks and the
  !spmd  taskid of the current instance.
  !
  ! Now read in a new value for N.  When it is 0, then you should depart.
  !
  read(20,*) n
  print *, "Number of approximation intervals = ", n
  !
  do while (n > 0)
    w = 1.0 / n
    sum = 0.0
    do i = 1,n
      sum = sum + f((i-0.5)*w)
    end do
    sum = sum * w
    err = sum - pi
    print *, "sum = ", sum, " err =", err
    !
    read (20,*) n
    print *, "Number of approximation intervals = ", n
  !
  end do
  !
  close (unit = 20)
  
contains
  real(DOUBLE) function f(x)
    implicit none
    real(DOUBLE), intent(in) :: x

    f = 4.0 / (1.0+x*x)
  end function f

end program karp
