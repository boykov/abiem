module modinteg
  use params
  use modint
  use modphi, only : phi
contains  
  
  include 'set_params.f90'

  subroutine inomp(i,nt)
    integer, intent(in) :: i, nt
    call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
    intphi_over(i) = folding(i,f,nt)
  end subroutine inomp

  double precision function f(x,i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f = phi(x,i,h**2) + 0
  end function f

  subroutine calcomp()
    use omp_lib
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call inomp(i,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp
  
end module modinteg
