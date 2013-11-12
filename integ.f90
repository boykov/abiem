module modinteg
  use params
  use modphi, only : phi
contains  
  
  include 'set_params.f90'

  double complex function approximateu(x)
    use omp_lib
    use dbsym
    double precision, intent(in), dimension(:) :: x
    integer :: i
    double complex, dimension(:), allocatable :: s

    allocate(s(numnodes))

    call OMP_SET_NUM_THREADS(4)


    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(v)
    do i=1,numnodes
       s(i) = intphi_over(i)*q(i)*Amn(x,node_coordinates(i,:),k)
    end do
    !$OMP END PARALLEL DO

    approximateu = sum(s(:))

    deallocate(s)
  end function approximateu

  double complex function matrixA(i,j)
    use omp_lib
    use dbsym
    integer, intent(in) :: i,j
    double precision :: sigm
    double precision, dimension(nd) :: x,y

    sigm = sigmaij(i,j)
    x = node_coordinates(i,:)
    y = node_coordinates(j,:)

    if (i .eq. j) then
       matrixA = (intphi_over(i)**2)*limA(sigm,k)
    else
       matrixA = intphi_over(i)*intphi_over(j)*Amn(x,y,k)
    end if
  end function matrixA

  double precision function sigmaij(i,j)
    integer, intent(in) :: i,j
    sigmaij = dsqrt(sigma(i)**2 + sigma(j)**2)
  end function sigmaij

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
  
  subroutine integrate(n,z,ip,nt)
    use dbsym
    integer, intent(in) :: ip,nt
    double precision, intent(in), dimension(:) :: n
    double precision, intent(in), dimension(:) :: z

    double precision rh,ph

    double precision, dimension(nd) :: x, y
    double precision, dimension(nd,nd) :: bt

    integer m,l,k,i,nthread

    integer Nk, iz, ik
    double precision jac
    nthread = nt + 1
    Nk = 4*Nz

    k = MINLOC(n,dim=1)

    do m=1,nd
       do l=1,nd
          bt(m,l) = beta(m,k,l,n,axes)
       end do
    end do

    do iz=1,Nz
       rh = centres(iz)

       do ik=1,Nk
          ph = (2.D0*PI/Nk)*ik
          do i=1,3
             x(i) = fx(i,bt,axes,h2,rh,ph,z)
             nodes(nthread,iz,ik,i) = x(i)
          end do

          jac = dsqrt(fjacobian(bt,axes,h2,rh,ph,z))
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*C(iz) * jac

       end do
    end do
    return
  end subroutine integrate

  double precision function folding(ip,f,nt)
    integer, intent(in) :: ip,nt
    interface
       function f(x,i)
         integer, intent(in) :: i
         double precision, intent(in), dimension(:) :: x
         double precision :: f
       end function f
    end interface

    integer Nk, iz, ik,nthread
    double precision gtmp, tmp

    Nk = 4*Nz
    nthread = nt + 1

    gtmp = 0
    do iz=1,Nz
       do ik=1,Nk
          tmp = f(nodes(nthread,iz,ik,:),ip)
          gtmp = gtmp + real(jacobian(nthread,iz,ik))*tmp
       end do
    end do
    folding = gtmp

  end function folding

  function getnode(nt,iz,ik)
    integer, intent(in) :: nt,iz,ik
    double precision, dimension(3) :: getnode ! TODO
    getnode(:) = nodes(nt,iz,ik,:)
  end function getnode

  double complex function foldingarr(ip,farr,nt)
    integer, intent(in) :: ip,nt
    double complex, dimension(:,:,:) :: farr
    integer Nk, iz, ik,nthread
    double complex gtmp, tmp

    Nk = 4*Nz
    nthread = nt + 1

    gtmp = (0,0)
    do iz=1,Nz
       do ik=1,Nk
          tmp = farr(nthread,iz,ik)
          gtmp = gtmp + real(jacobian(nthread,iz,ik))*tmp
       end do
    end do
    foldingarr = gtmp

  end function foldingarr

end module modinteg
