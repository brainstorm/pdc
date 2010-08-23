program prog1
  use omp_smeds
  implicit none

  integer :: i,j,thrd_id,old_maxtrds
  real, allocatable:: numbers(:)
  logical:: old_dynamic
  double precision:: val1, sum

  integer :: org_number_of_threads

!!! This set of exercises aims at showing some basic OpenMP directives and how
!!! they affect the order (and correctness) of execution

  i=0
  print *, "This is the first line in the program"
  !$ i=1
  if(i==0)then
    print *, "Oooops!  Something is wrong. This program is not OpenMP enabled."
    print *, "         Perhaps you have forgotten some compiler switches???"
    stop
  end if

  allocate(numbers(omp_get_max_threads()))
  numbers(:)=0
  print *,'size=',size(numbers)
!!! ------------------------------------------------------------

  write (*,"(//a/)") "   ----- Exercise 1 -------"
!!! This exercise tries to illustrate the simplest parallel OpenMP
!!! program that can be written. Run it several times. At some occasions
!!! the printed output is "wrong". Why is that? Correct the program so
!!! that it is "correct" 

!$omp parallel private(thrd_id) 
  print *, "This is a line in the parallel section"
  !$  thrd_id=omp_get_thread_num()
  call random_number(numbers(thrd_id+1))
  print *, "My value of thrd_id is",thrd_id
!$omp end parallel

#if EXAMPLE == 1
  write(*,"(/a/)") "   ----- End of exercise 1 ------"
  stop 
#endif

  write (*,"(//a/)") "   ----- Exercise 2 -------"
!!! Use different values for the environment variable OMP_NUM_THREADS and
!!! run this program. Now, modify the program so that this section always 
!!! runs on as many threads as there are processors. Also make sure that
!!! the num_threads is set back to its original value afterwards
  ! ----> Insert code here
  !$  org_number_of_threads = omp_get_max_threads()
  !$ call omp_set_num_threads(omp_get_num_procs())
!$omp parallel 
  print *,"I am thread",omp_get_thread_num()," of a total of",omp_get_num_threads()
!$omp end parallel
  ! ----> Insert code here too
  !$ call omp_set_num_threads(org_number_of_threads)
#if EXAMPLE == 2
  write(*,"(/a/)") "   ----- End of exercise 2 ------"
  stop
#endif

!!! ------------------------------------------------------------

  write (*,"(//a/)") "   ----- Exercise 3 -------"
!!! Run this example several times. You will notice that the master output 
!!! occasionally comes before and occasionally after the single thread
!!! statement. Change the program so that the master task always 
!!! writes its statement before the single task does
!$omp parallel private(thrd_id)
  thrd_id=omp_get_thread_num()

!$omp sections
  print *,"This is task A running in thread",thrd_id
!$omp section
  print *,"This is task B running in thread",thrd_id
!$omp section
  print *,"This is task C running in thread",thrd_id
!$omp end sections nowait

!$omp master
  print *,"This is the master thread running at thread",thrd_id
!$omp end master
!$omp barrier
!$omp single
  print *,"This task was executed at thread",thrd_id
!$omp end single nowait

!$omp end parallel

#if EXAMPLE == 3
  write(*,"(/a/)") "   ----- End of exercise 3 ------"
  stop
#endif

!!! ------------------------------------------------------------

  write (*,"(//a/)") "   ----- Exercise 4 -------"

!!! schedule dynamic example
!!! When scheduling is set to dynamic the operating system is free to select
!!! the number of threads to use in a parallel section. It can not use more
!!! than omp_get_num_threads(), but it is allowed to use fewer.

!!! A system that allows for dynamic scheduling and that has a defensive approach
!!! of issuing the number of processes might show interesting output in this
!!! example. The current (August, 2000) installation does not however.

  old_dynamic=omp_get_dynamic()
  call omp_set_dynamic(.true.)
!$omp parallel
  print *,"I am thread",omp_get_thread_num()," of a total of",omp_get_num_threads()
!$omp single
  print *,"Report from thread",omp_get_thread_num()
  print *,"The maximum allowed number of threads is",omp_get_max_threads()
  print *,"In adddition, the dynamic scheduling evaluates to",omp_get_dynamic()
!$omp end single
  
!$omp end parallel

  call omp_set_dynamic(old_dynamic)

#if EXAMPLE == 4
  write(*,"(/a/)") "   ----- End of exercise 4 ------"
  stop
#endif

!!! ------------------------------------------------------------
  write (*,"(//a/)") "   ----- END -------"

end program prog1
