!!! This program generates some vectors and matrices
!!! It manipulates them and finally computes some global
!!! things that are printed


program Prog2
  use omp_smeds
  implicit none
  integer, parameter:: MATSIZE=900
  double precision, dimension(:,:), allocatable :: Mat_A, Mat_B, Mat_C, Mat_D
  double precision :: x,y,scal,Sum_A,Sum_B,Sum_C,Sum_D, xA, xB, xC
  common scal,x,y
  integer :: i,j,k

  real(8),external :: rtc
  real(8) :: dtime,dtime2
  external dummy

  allocate(Mat_A(MATSIZE,MATSIZE),Mat_B(MATSIZE,MATSIZE),Mat_C(MATSIZE,MATSIZE))
  ! Generate three non-random matrices

  write(*,"(//a/)") "   ----- Exercise 1 ------"
!!! The code below generates three matrices. Try to think of a way in which
!!! this can be made parallel in any way. Make sure that the printed output
!!! x is correct in your parallel version

  dtime=rtc()

!$omp parallel private(xA, xB, xC, i,j) shared(x, Mat_A, Mat_B, Mat_C)
!$omp sections
  xA=0.35d0
  do j=1,MATSIZE
    do i=1,MATSIZE
      xA=1-fraction(sin(xA))
      Mat_A(i,j)=xA
    end do
  end do

!$omp section
  xB=0.68d0
  do j=1,MATSIZE
    do i=1,MATSIZE
      xB=1-fraction(sin(xB))
      Mat_B(i,j)=xB
    end do
  end do
      
!$omp section
  xC=0.24d0
  do j=1,MATSIZE
    do i=1,MATSIZE
      xC=1-fraction(sin(xC))
      Mat_C(i,j)=xC
    end do
  end do
  x = xC
!$omp end sections
!$omp end parallel

  dtime=rtc()-dtime
  Sum_A=check_sum(Mat_A)
  Sum_B=check_sum(Mat_B)
  Sum_C=check_sum(Mat_C)

  print *," The check sum of the matrices evaluates to:"
  print 100,"A",Sum_A
  print 100,"B",Sum_B
  print 100,"C",Sum_C
  print 110,dtime*1d3

  print 101,"The variable x evaluates to",x
100 format("Sum of matrix ",a,g25.16)
101 format(a,g37.25)
110 format("Time for the exercise: ",f9.1,"ms")
#if EXAMPLE == 1
  write(*,"(/a/)") "   ----- End of exercise 1 ------"
  stop 
#endif


  write(*,"(//a/)") "   ----- Exercise 2 ------"
!!! This code makes a simple attempt at computing a matrix multiply. Try
!!! to parallelize it without changing the results (more than negligible)
!!! In this exercise parallelize the outer-most loop

  allocate(Mat_D(MATSIZE,MATSIZE))
  dtime=rtc()
 
!$omp parallel do private(i,j,k) shared(Mat_A, Mat_B, Mat_D) lastprivate(scal)
  do i=1,MATSIZE
    do j=1,MATSIZE
      scal=0.0d0
      do k=1,MATSIZE
        scal=scal+Mat_A(i,k)*Mat_B(k,j)
      end do
      Mat_D(i,j)=scal
    end do
  end do
  dtime=rtc()-dtime

  Sum_D=check_sum(Mat_D)
  print *," The check sum of the matrices evaluates to:"
  print 100,"D",Sum_D
  print 101,"The value of scal is:",scal
  print 110,dtime*1d3

#if EXAMPLE == 2
  write(*,"(/a/)") "   ----- End of exercise 2 ------"
  stop 
#endif

  write(*,"(//a/)") "   ----- Exercise 3 ------"
!!! This code makes a simple attempt at computing a matrix multiply. Try
!!! to parallelize it without changing the results (more than negligible)
!!! In this exercise parallelize the second outer-most loop

  dtime=rtc()
  
  do i=1,MATSIZE
!$omp parallel do private(j,k) shared(Mat_A, Mat_B, Mat_D) lastprivate(scal)
    do j=1,MATSIZE
      scal=0.0d0
      do k=1,MATSIZE
        scal=scal+Mat_A(i,k)*Mat_B(k,j)
      end do
      Mat_D(i,j)=scal
    end do
  end do
  dtime=rtc()-dtime

  Sum_D=check_sum(Mat_D)
  print *," The check sum of the matrices evaluates to:"
  print 100,"D",Sum_D
  print 101,"The value of scal is:",scal
  print 110,dtime*1d3

#if EXAMPLE == 3
  write(*,"(/a/)") "   ----- End of exercise 3 ------"
  stop 
#endif


  write(*,"(//a/)") "   ----- Exercise 4 ------"
!!! This code shows a much better way of doing the matrix multiply
!!! What is the difference between the serial and threaded programs?

  dtime=rtc()
  call DGEMM('N','N',MATSIZE,MATSIZE,MATSIZE,1.0d0,Mat_A,MATSIZE,&
       &  Mat_B,MATSIZE,0.0d0,Mat_D,MATSIZE)
  dtime=rtc()-dtime

  Sum_D=check_sum(Mat_D)
  print *," The check sum of the matrices evaluates to:"
  print 100,"D",Sum_D
  print 110,dtime*1d3

#if EXAMPLE == 4
  write(*,"(/a/)") "   ----- End of exercise 4 ------"
  stop 
#endif

  write(*,"(//a/)") "   ----- Exercise 5 ------"
!!! In this example we will sum the element in the matrix. 
!!! Try to find a way to parallelize the summation made in
!!! the dtime2 timing loop
  print *,"First a dummy round to warm up the cache"
  print *,"This first summation might take longer time than the next one."
  call dummy(Sum_C)
  call dummy(Mat_C)
  dtime=rtc()
  Sum_C=sum(Mat_C)
  dtime=rtc()-dtime
  call dummy(Sum_C)
  call dummy(Mat_C)
  print 110,dtime*1d3

  dtime=rtc()
  Sum_C=sum(Mat_C)
  dtime=rtc()-dtime

  !  --->  This is the section of the code to parallelize
  dtime2=rtc()
  scal=0.0d0
  do j=1,MATSIZE
    do i=1,MATSIZE
      scal=scal+Mat_C(i,j)
    end do
  end do
  dtime2=rtc()-dtime2

  print *," The check sum of the matrices evaluates to:"
  print 100,"C (version1)",Sum_C
  print 100,"C (version2)",scal
  print 100,"The relative error is (binary units)", &
       & abs(Sum_C-scal)/(Sum_C-nearest(Sum_C,1.0))
  print 110,dtime*1d3
  print 110,dtime2*1d3

#if EXAMPLE == 5
  write(*,"(/a/)") "   ----- End of exercise 5 ------"
  stop 
#endif


contains
  function check_sum(Mat)
    implicit none
    double precision :: check_sum
    double precision :: Mat(:,:)
    check_sum=sum(Mat)
  end function check_sum
end program Prog2
