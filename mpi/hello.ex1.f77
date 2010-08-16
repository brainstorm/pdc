C Program hello.ex1.f
C Parallel version using MPI calls
C Modified from basic version so that workers send back
C a message to the master, who prints out a message
C for each worker
	program hello
	include 'mpif.h'
	integer rank,size,mpierr,tag,status(MPI_STATUS_SIZE)
	integer wrank
	character(12) message, inmsg

	call MPI_INIT(mpierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierr)
        tag = 100

	if(rank .eq. 0) then
	  message = 'Hello, world'
	  do i=1,size-1
             call MPI_SEND(message,12,MPI_CHARACTER,i,tag,
     .                     MPI_COMM_WORLD,mpierr)
          enddo
          write(6,*)'process',rank,':',message
          do i=1,size-1
             call MPI_RECV(wrank,1,MPI_INTEGER,MPI_ANY_SOURCE,tag,
     .                     MPI_COMM_WORLD,status,mpierr)
             write(6,*) 'process', wrank, ':Hello, back'
          enddo
       else
	  call MPI_RECV(inmsg,12,MPI_CHARACTER,0,tag,
     .                  MPI_COMM_WORLD,status,mpierr)
          call MPI_SEND(rank,1,MPI_INTEGER,0,tag,
     .                  MPI_COMM_WORLD,mpierr)
       endif
       call MPI_FINALIZE(mpierr)
       end 

