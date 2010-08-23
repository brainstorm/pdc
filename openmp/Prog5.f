module realtypes
  integer, parameter:: SINGLE=kind(1.0), DOUBLE=kind(1.0d0)
end module realtypes

program prog5
  use realtypes
  !$ use omp_smeds
  implicit none

  integer, parameter:: N=2123123,NBINS=16
  real(DOUBLE), parameter :: PI=3.141592653589793238462643383279506d0
  integer :: i,j,dN,tid
  real(DOUBLE), allocatable, dimension(:,:) :: A
  real(DOUBLE), allocatable, dimension(:) :: blocksum
  real(DOUBLE), dimension(NBINS) :: binval
  integer, dimension(NBINS) :: bins
  real(DOUBLE) :: dw
  integer :: Nthrds,start,end,MaxNthrds

  real(8) :: dtime, rtc
  external :: rtc

  MaxNthrds=1
  !$ MaxNthrds=omp_get_max_threads()

  allocate(A(N,2),blocksum(0:MaxNthrds-1))

  dw=PI/real(NBINS,kind=DOUBLE)
  do i=1,NBINS-1   !!! Set up the upper limit for each bin
    binval(i)=i*dw
  end do
  binval(NBINS)=huge(binval(NBINS))
  print *,"The binvals: ",binval

  dw=PI/real(N,kind=DOUBLE)

  do i=1,N  ! Initiate the array with values in a "random" order
    A(i,1)=(1+modulo(i*1077,N))*dw
  end do

  print *

!!!-------------------------------------------------------

  bins(:)=0
  dtime=rtc()
  A_loop: do i=1,N
    do j=1,NBINS
      if(A(i,1)<binval(j))then 
        bins(j)=bins(j)+1 
        cycle A_loop
      endif
    end do
  end do A_loop 
  dtime=rtc()-dtime
  print *,"First summation loop took",dtime*1d3,"ms"
  print *,"The bins are:",bins
  print * 

!!!-------------------------------------------------------

  bins(:)=0
  dtime=rtc() 
  A2_loop: do i=1,N
    do j=1,NBINS
      if(A(i,1)<binval(j))then 
        bins(j)=bins(j)+1 
        cycle A2_loop
      endif
    end do
  end do A2_loop 
  dtime=rtc()-dtime
  print *,"Second summation loop took",dtime*1d3,"ms"
  print *,"The bins are:",bins
  print *

!!!-------------------------------------------------------

  bins(:)=0
  dtime=rtc() 
  Aloop: do i=1,N
    do j=1,NBINS
      if(A(i,1)<binval(j))then 
        bins(j)=bins(j)+1
        cycle Aloop
      endif
    end do
  end do Aloop 
  dtime=rtc()-dtime
  print *,"Third summation loop took",dtime*1d3,"ms"
  print *,"The bins are:",bins
  print *

!!!-------------------------------------------------------

  bins(:)=0
  dtime=rtc() 
  AA_loop: do i=1,N
    do j=1,NBINS
      if(A(i,1)<binval(j))then
        bins(j)=bins(j)+1
        cycle AA_loop
      endif
    end do
  end do AA_loop 
  dtime=rtc()-dtime
  print *,"Fourth summation loop took",dtime*1d3,"ms"
  print *,"The bins are:",bins
  print *

!!!-------------------------------------------------------


end program prog5
