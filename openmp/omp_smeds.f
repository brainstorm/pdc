!!!$----------------------------------------------------
!!!$
!!!$   This is a suggested module to declare the OpenMP
!!!$   functions and subroutine interfaces. This 
!!!$   interface also declares a type for the OpenMP
!!!$   lock that hopefully is portable between systems.
!!!$   
!!!$   A lock is in a user program is declared as:
!!!$      use omp_smeds
!!!$      type(omp_lock_t):: lock
!!!$
!!!$  The declaration is such that it is possible to 
!!!$  construct allocatable arrays of locks
!!!$
!!!$      type(omp_lock_t), allocatable:: lock_array(:)
!!!$
!!!$  Everything is provided as is and with no guarantees
!!!$  Nils Smeds <smeds@pdc.kth.se>, 2000-08-28
!!!$
!!!$----------------------------------------------------

module omp_smeds_types
  type, public :: omp_lock_t
    sequence   !!! Based on the standards POINTER(IL,LOCK)
    integer, pointer :: omp_lock_t_ptr
  end type
end module omp_smeds_types

module omp_smeds
  use omp_smeds_types
  interface
    subroutine omp_set_num_threads(NP)
      integer:: NP
    end subroutine omp_set_num_threads

    function omp_get_num_threads()
      integer:: omp_get_num_threads
    end function omp_get_num_threads

    function omp_get_max_threads()
      integer:: omp_get_max_threads
    end function omp_get_max_threads

    function omp_get_thread_num()
      integer:: omp_get_thread_num
    end function omp_get_thread_num

    function omp_get_num_procs()
      integer:: omp_get_num_procs
    end function omp_get_num_procs

    subroutine omp_set_dynamic(flag)
      logical:: flag
    end subroutine omp_set_dynamic

    function omp_get_dynamic()
      logical:: omp_get_dynamic
    end function omp_get_dynamic

    function omp_in_parallel()
      logical:: omp_in_parallel
    end function omp_in_parallel

    subroutine omp_set_nested(flag)
      logical:: flag
    end subroutine omp_set_nested

    function omp_get_nested()
      logical:: omp_get_nested
    end function omp_get_nested

    subroutine omp_init_lock(lock)
      use omp_smeds_types
      type(omp_lock_t) :: lock
    end subroutine omp_init_lock

    subroutine omp_destroy_lock(lock)
      use omp_smeds_types
      type(omp_lock_t) :: lock
    end subroutine omp_destroy_lock

    subroutine omp_unset_lock(lock)
      use omp_smeds_types
      type(omp_lock_t) :: lock
    end subroutine omp_unset_lock

    function omp_test_lock(lock)
      use omp_smeds_types
      type(omp_lock_t) :: lock
      logical :: omp_test_lock
    end function omp_test_lock

  end interface

end module omp_smeds
