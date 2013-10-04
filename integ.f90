module modinteg
  use params
  use intphi, only : inomp
contains  
  
  include 'set_params.f90'

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
