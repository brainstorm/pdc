! Program hello.ex2.f
!
! Parallel version using MPI calls.
! Modified from basic version so that workers send back a message to the 
! master, who prints out a message for each worker.  In addition, the 
! master now sends out two messages to each worker, with two different 
! tags, and the worker receives the messages in reverse order.
!
! Note that this solution works only because the messages are small,
! and can fit into buffers.  A later talk will provide details on
! how buffers are used in MPI_SEND and MPI_RECEIVE,
!
program hello
  implicit none
  integer, parameter:: DOUBLE=kind(1.0d0), SINGLE=kind(1.0)
  include "mpif.h"
  character(LEN=12) :: inmsg,message
  integer :: i,mpierr,rank,size,tag,tag2,wrank
  integer, dimension(MPI_STATUS_SIZE) :: status
  !
  call MPI_Init(mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD,size,mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierr)
  tag = 100
  tag2 = 200
  !
  if (rank == 0) then
    message = "Hello, world"
    do i = 1,size-1
      call MPI_Send(message,12,MPI_CHARACTER,i,tag,MPI_COMM_WORLD,mpierr)
      call MPI_Send(message,12,MPI_CHARACTER,i,tag2,MPI_COMM_WORLD,mpierr)
    end do
    write(*,*) "process", rank, ":", message
    do i = 1,size-1
      call MPI_Recv(wrank,1,MPI_INTEGER,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD, &
                    status,mpierr)
      write(*,*) "process", wrank, ":Hello, back"
    end do
  else
    call MPI_Recv(inmsg,12,MPI_CHARACTER,0,tag2,MPI_COMM_WORLD,status,mpierr)
    call MPI_Recv(inmsg,12,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,status,mpierr)
    call MPI_Send(rank,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,mpierr)
  end if
  call MPI_Finalize(mpierr)
end program hello
!
