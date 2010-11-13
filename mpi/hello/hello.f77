C Program hello.f
	program hello
	include 'mpif.h'
	integer rank,size,mpierr,tag,status(MPI_STATUS_SIZE)
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
       else
	  call MPI_RECV(inmsg,12,MPI_CHARACTER,0,tag,
     .                  MPI_COMM_WORLD,status,mpierr)
	  write(6,*)'process',rank,':',inmsg
       endif
       call MPI_FINALIZE(mpierr)
       end 

