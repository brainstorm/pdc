!*******************************************************************    
!   twod.f - a solution to the Poisson problem using Jacobi             
!   interation on a 1-d decomposition                                   
!                                                                       
!   The size of the domain is read by processor 0 and broadcast to      
!   all other processors.  The Jacobi iteration is run until the        
!   change in successive elements is small or a maximum number of       
!   iterations is reached.  The difference is printed out at each       
!   step.                                                               
!
!*******************************************************************    
    Program main
      Implicit none
!                                                                       !
! ..
! .. Local Scalars ..
      Real (8) :: diffnorm=1.0, rtime, tolerance = 1.E-06
      Integer :: it=0, nx, ny, t1, t2, count_rate
! ..
! .. Local Arrays ..
      Real (8), Allocatable, Dimension (:,:) :: a, b, f
! ..
! .. External Functions ..
!      Real (8), External :: diff
!IBM      Integer, External :: mclock
!IBM      Integer:: c1, ctime, 
      Real:: c1, ctime, ctime_s
! ..

!                                                                       
!         Get the size of the problem                                   
!                                                                       
      print *, 'Enter nx'                                          
      read *, nx                                                   
      ny = nx

!
! Allocate node local arrays
!
      Allocate(a (0:nx+1,0:ny+1), b (0:nx+1,0:ny+1), f (0:nx+1,0:ny+1))
!                                                                       
! Initialize the right-hand-side (f) and the initial solution guess (a) 
!                                                                       
      Call twodinit(a,b,f,nx)

!                                                                       
! Actually do the computation.  Note the use of a collective operation t
! check for convergence, and a do-loop to bound the number of iterations
!                                                                       
      Call system_clock(t1, count_rate)
!IBM      c1 = mclock()
      call CPU_TIME(c1)
      Do While ( diffnorm > tolerance)
        it = it + 1
        Call sweep2d(a, f, nx, b)
        Call sweep2d(b, f, nx, a)

        diffnorm = diff(a(1:nx,1:ny), b(1:nx,1:ny))

        If (2*it<=1000 .and. mod(2*it,100)==0) Write(*,*) 2*it,&
             & ' Difference is ', diffnorm
        If (2*it>1000 .and. mod(2*it,1000)==0) Write(*,*) 2*it,&
             & ' Difference is ', diffnorm
      End Do

!IBM      ctime = mclock() - c1
!IBM      ctime_s = ctime / 100.d0
      call CPU_TIME(ctime)
      ctime_s = ctime - c1

      Call system_clock(t2)
      rtime = dble(t2 - t1)/dble(count_rate)
      Write(*,*)
      Write(*,*) 'Converged after ', 2*it, ' Iterations'
      Write(*,*) 2*it, ' Difference is ', diffnorm
      Write(*,'(a, f12.2)') ' Cpu time (s):  ', ctime_s
      Write(*,'(a, f12.2)') ' Real time (s): ', rtime

!.. Write solution to  file twod.out        
      open(unit=20,file='twod.out',status='UNKNOWN')
      write(20,*) b(0:nx+1,0:ny+1)
!                            

    contains

!                                                                       
!  The rest of the 1-d program                                          
!                                                                       
    Function diff(a,b) 
!                                                                       
! .. Function Return Value ..
      Double Precision :: diff
! ..
! .. Array Arguments ..
      Real (8), Dimension (:,:) :: a, b
! ..
! .. Local scalars
      integer :: i,j
      real(8) :: diffsum, tmp

      diffsum=0.0d0
     !$omp parallel do default(none) private(j, i, tmp) shared(a, b) reduction(+: diffsum)
      do j=1,ubound(a,DIM=2)
        do i=1,ubound(a,DIM=1)
          tmp=a(i,j)-b(i,j)
          diffsum = diffsum + tmp*tmp
        end do
      end do

      diff=diffsum
      Return
    End Function diff

    Subroutine twodinit(a, b, f, nx)
!                                                                       
! .. Scalar Arguments ..
      Integer :: nx, ny, i,j
! ..
! .. Array Arguments ..
      Real (8), Dimension (0:nx+1,0:nx+1) :: a, b, f

! ..
!                                                                       
      ny = nx

      a = 0.D0
      b = 0.D0
      f = 0.D0
!                                                                       
!    Handle boundary conditions                                         
!                                                                       
      a(0:0,0:ny+1) = 1.D0
      b(0:0,0:ny+1) = 1.D0

      a(0:nx+1,0:0) = 1.D0
      b(0:nx+1,0:0) = 1.D0

!                                                                       
      Return
    End Subroutine twodinit

!                                                                       
! Perform a Jacobi sweep for a 1-d decomposition.                       
! Sweep from a into b                                                   
!                                                                       
    Subroutine sweep2d(a, f, nx, b)
!                                                                       
! .. Scalar Arguments ..
      Integer :: nx, ny
! ..
! .. Array Arguments ..
      Real (8), Dimension (0:nx+1,0:nx+1) :: a, b, f

! ..
! .. Local Scalars ..
      Real (8) :: h
      Integer :: i, j
! ..
! .. Intrinsic Functions ..
      Intrinsic dble
! ..
!                    
      ny = nx
      h = 1.0d0/dble(nx+1)
      
      !$omp parallel do default(none) private(j, i) shared(a, b, nx, ny, h, f)
      do j=1,ny
        do i=1,nx
          b(i,j) = 0.25*(a(i-1,j)+a(i,j+1)+a(i,j-1)+a(i+1,j)) - h*h*f(i,j)
        End do
      End do

      Return
    End Subroutine sweep2d

                                           
    End Program main





