    PROGRAM example1
      IMPLICIT NONE

      INTEGER, PARAMETER :: maxn = 14000
      INTEGER :: i
      REAL (8) :: a, b, c

      a = 0
      DO i = 1, maxn
        call random_number(b)
        call random_number(c)
        a = a + b + c
      END DO

      PRINT *, 'Number of iterations:', maxn, '  Result: ', a

    END PROGRAM example1
