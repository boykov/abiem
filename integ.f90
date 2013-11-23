module modinteg
  use params
  use modphi, only : phi
contains  
  
  include 'set_params.f90'

  double precision function test_fast()
    use fast_dbsym
    test_fast = foo2()
  end function test_fast

  double complex function matrixA2(i,j)
    use dbsym
    integer, intent(in) :: i,j
    integer :: l
    double precision :: sigm, nstar, s
    double precision, dimension(3) :: x,y

    nstar = sum(normal_coordinates(i,:)*(node_coordinates(i,:)-node_coordinates(j,:))*intphi_over(j))

    sigm = sigmaij(i,j)
    x = node_coordinates(i,:)
    y = node_coordinates(j,:)

    if (i .eq. j) then
       matrixA2 = (gauss(j))*(intphi_over(i))
    else
       ! print *, Bmn(x,y,k)
       matrixA2 = intphi_over(i)*nstar*Bmn(x,y,k)
    end if
  end function matrixA2

  subroutine setgauss()
    use dbsym
    integer :: l,j
    double precision :: sigm, nstar, s
    double precision, dimension(3) :: x,y

    do j=1,numnodes
       s = 0.0
       s2 = 0.0
       do l=1,numnodes
          if (l .ne. j) then
             nstar = sum(normal_coordinates(l,:)*(node_coordinates(l,:)-node_coordinates(j,:))*intphi_over(l))
             sigm = sigmaij(l,j)
             x = node_coordinates(l,:)
             y = node_coordinates(j,:)
             s = s + nstar/(4*PI*norm(node_coordinates(l,:)-node_coordinates(j,:))**3)
          end if
       end do
       gauss(j) = s
    end do
  end subroutine setgauss

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

  double complex function approximateu_sigm(x)
    use omp_lib
    use dbsym
    double precision, intent(in), dimension(:) :: x
    double precision :: sigm
    integer :: i
    double complex, dimension(:), allocatable :: s

    allocate(s(numnodes))

    call OMP_SET_NUM_THREADS(4)


    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(v)
    do i=1,numnodes
       sigm = sigma(i)
       s(i) = intphi_over(i)*q(i)*A(x,node_coordinates(i,:),sigm,k)
    end do
    !$OMP END PARALLEL DO

    approximateu_sigm = sum(s(:))

    deallocate(s)
  end function approximateu_sigm

  double complex function matrixA_sigm(i,j)
    use omp_lib
    use dbsym

    integer, intent(in) :: i,j
    double precision :: sigm
    double precision, dimension(nd) :: x,y

    sigm = sigmaij(i,j)
    x = node_coordinates(i,:)
    y = node_coordinates(j,:)

    if (i .eq. j) then
       matrixA_sigm = (intphi_over(i)**2)*limA(sigm,k)
    else
       matrixA_sigm = intphi_over(i)*intphi_over(j)*A(x,y,sigm,k)
    end if
  end function matrixA_sigm

  double complex function vectorb(i)
    use dbsym
    integer, intent(in) :: i

    vectorb = cdexp((0,1)*k*node_coordinates(i,3))*intphi_over(i)

  end function vectorb

  double complex function vectorb2(i)
    use dbsym
    integer, intent(in) :: i

    vectorb2 = (0,1)*k*normal_coordinates(i,3)*cdexp((0,1)*k*node_coordinates(i,3))*intphi_over(i)

  end function vectorb2

  double precision function sigmaij(i,j)
    integer, intent(in) :: i,j
    sigmaij = dsqrt(sigma(i)**2 + sigma(j)**2)
  end function sigmaij

  subroutine inomp(i,nt)
    use dbsym
    integer, intent(in) :: i, nt

    ptr_jacobian => fjacobian
    call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
    intphi_over(i) = folding(i,f,nt)
  end subroutine inomp

  subroutine inomp2(i,nt)
    use dbsym
    integer, intent(in) :: i, nt
    double complex :: zero

    zero = 0
    ptr_jacobian => fjacobian2
    call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
    intphi_under(i) = folding(i,f,nt)
  end subroutine inomp2

  double precision function f2(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f2 = phi(x,i,h**2)/norm(x)
  end function f2

  double precision function f(x,i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f = phi(x,i,h**2) + 0
  end function f

  subroutine calcomp2()
    use omp_lib
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call inomp2(i,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp2

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

          jac = dsqrt(ptr_jacobian(bt,axes,h2,rh,ph,z))
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*C(iz) * jac

       end do
    end do
    return
  end subroutine integrate

  double precision function folding(ip,f,nt)
    use dbsym
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
          gtmp = gtmp + realpart(jacobian(nthread,iz,ik))*tmp
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
