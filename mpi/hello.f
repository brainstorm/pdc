! Program hello.f
program hello
  implicit none
  integer, parameter:: DOUBLE=kind(1.0d0), SINGLE=kind(1.0)
  include "mpif.h"
  character(LEN=12) :: inmsg,message
  integer :: i,mpierr,rank,size,tag
  integer, dimension(MPI_STATUS_SIZE) :: status
  !
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD,size,mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierr)
  tag = 100
  !
  if (rank == 0) then
    message = "Hello, world"
    do i = 1,size-1
      call MPI_Send(message,12,MPI_CHARACTER,i,tag,MPI_COMM_WORLD,mpierr)
    end do
    write(*,*) "process", rank, ":", message
  else
    call MPI_Recv(inmsg,12,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,status,mpierr)
    write(*,*) "process", rank, ":", inmsg
  end if
  call MPI_Finalize(mpierr)
end program hello
!
