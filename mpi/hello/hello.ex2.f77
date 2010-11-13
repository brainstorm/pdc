C Program hello.ex2.f
C
C Parallel version using MPI calls.
C Modified from basic version so that workers send back a message to the 
C master, who prints out a message for each worker.  In addition, the 
C master now sends out two messages to each worker, with two different 
C tags, and the worker receives the messages in reverse order.
C
C Note that this solution works only because the messages are small,
C and can fit into buffers.  A later talk will provide details on
C how buffers are used in MPI_SEND and MPI_RECEIVE,
C
	program hello
	include 'mpif.h'
	integer rank,size,mpierr,tag,status(MPI_STATUS_SIZE)
	integer wrank, tag2
	character(12) message, inmsg

	call MPI_INIT(mpierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierr)
        tag = 100
        tag2 = 200

	if(rank .eq. 0) then
	  message = 'Hello, world'
	  do i=1,size-1
             call MPI_SEND(message,12,MPI_CHARACTER,i,tag,
     .                     MPI_COMM_WORLD,mpierr)
             call MPI_SEND(message,12,MPI_CHARACTER,i,tag2,
     .                     MPI_COMM_WORLD,mpierr)
          enddo
          write(6,*)'process',rank,':',message
          do i=1,size-1
             call MPI_RECV(wrank,1,MPI_INTEGER,MPI_ANY_SOURCE,tag,
     .                     MPI_COMM_WORLD,status,mpierr)
             write(6,*) 'process', wrank, ':Hello, back'
          enddo
       else
	  call MPI_RECV(inmsg,12,MPI_CHARACTER,0,tag2,
     .                  MPI_COMM_WORLD,status,mpierr)
	  call MPI_RECV(inmsg,12,MPI_CHARACTER,0,tag,
     .                  MPI_COMM_WORLD,status,mpierr)
          call MPI_SEND(rank,1,MPI_INTEGER,0,tag,
     .                  MPI_COMM_WORLD,mpierr)
       endif
       call MPI_FINALIZE(mpierr)
       end 

