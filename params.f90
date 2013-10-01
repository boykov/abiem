module params
  ! number of dimensions
  integer, parameter :: nd = 3
  double precision, parameter :: PI = 3.14159265358979324D0
    integer :: Nz
  double precision, dimension(nd) :: axes
  double precision :: h2
  double precision, dimension(:), allocatable       :: centres
  double precision, dimension(:), allocatable       :: C
  double complex, dimension(:,:,:), allocatable   :: jacobian
  double precision, dimension(:,:,:,:), allocatable :: nodes
contains  
  double precision function asqrt(x)
    double precision, intent(in) :: x
    asqrt = dsqrt(abs(x))
  end function asqrt

  integer function dn(n,m)
    integer :: n,m
    dn = 0
    if (n .eq. m) dn = 1
  end function dn
  
  PURE SUBROUTINE Insertion_Sort(a)
    double precision, INTENT(in out), DIMENSION(:) :: a
    double precision :: temp
    INTEGER :: i, j

    DO i = 2, SIZE(a)
       j = i - 1
       temp = a(i)
       DO WHILE (j>=1 .AND. a(j)>temp)
          a(j+1) = a(j)
          j = j - 1
       END DO
       a(j+1) = temp
    END DO
  END SUBROUTINE Insertion_Sort
  
  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed
end module params
