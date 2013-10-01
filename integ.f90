module modinteg
  use params
  integer :: max_neighbors
  integer :: numnodes
  double precision :: h

  double precision, pointer :: node_coordinates(:,:)
  double precision, pointer :: normal_coordinates(:,:)
  double precision, pointer :: nstroke_coordinates(:,:)

  ! double precision, pointer :: axes(:)

  integer, pointer :: node_neighbors2(:,:)
  integer, pointer :: node_neighbors1(:,:)

  double precision, dimension(:), allocatable   :: intphi_over
  double precision :: area
  
contains
  double precision function calcarea()
    use omp_lib
    use modint, only : integrate, folding, foldingarr
    integer i, nt,iz,ik
    double precision tmp

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
       intphi_over(i)       = folding(i,f,nt)
    end do
    !$OMP END PARALLEL DO
    area = sum(intphi_over)
    calcarea = area
  end function calcarea

  double precision function f(x,i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(3) :: x
    f = 1 ! phi(x,i,h**2)
  end function f
  
end module modinteg
