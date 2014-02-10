module modinteg
  use params
  use modphi, only : phi, deltah
  integer :: i_tmp, j_tmp
  double precision :: y_tmp(3)
contains

  include 'set_params.f90'
  include 'kernels.f90'

  double precision function fAre(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    fAre = phi(x,i,hval**2)*realpart(Amn(y_tmp,x(:) + node_coordinates(i,:),k_wave))
  end function fAre

  double precision function fAim(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    fAim = phi(x,i,hval**2)*imagpart(Amn(y_tmp,x(:) + node_coordinates(i,:),k_wave))
  end function fAim

  double precision function f2(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f2 = phi(x,i,hval**2)/norm(x(:) + node_coordinates(i,:) - node_coordinates(j_tmp,:))
  end function f2

  double precision function f(x,i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f = phi(x,i,hval**2)
  end function f

  double precision function f3(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f3 = (1 - deltah(x(:)  + node_coordinates(i,:) - node_coordinates(j_tmp,:),hval))*phi(x,i,hval**2)/norm(x(:) + node_coordinates(i,:) - node_coordinates(j_tmp,:))
  end function f3

  double precision function f4(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f4 = (deltah(x(:) + node_coordinates(i,:) - node_coordinates(j_tmp,:),hval))*phi(x(:) + node_coordinates(i,:) - node_coordinates(i_tmp,:) ,i_tmp,hval**2)
  end function f4

  double precision function test_fast()
    use fast_dbsym
    test_fast = foo2()
  end function test_fast

!                 _
!        ___  ___| |_
!       / __|/ _ \ __|
!       \__ \  __/ |_
!       |___/\___|\__|
!

  subroutine setgauss()
    use dbsym
    integer :: l,j
    double precision :: sigm, nstar, s
    double precision, dimension(3) :: x,y
    double complex :: s2, s3

    do j=1,numnodes
       s = 0.0
       s2 = 0.0
       s3 = 0.0
       do l=1,numnodes
          if (l .ne. j) then
             nstar = sum(normal_coordinates(l,:)*(node_coordinates(l,:)-node_coordinates(j,:))*intphi_over(l))
             sigm = sigmaij(l,j)
             x = node_coordinates(l,:)
             y = node_coordinates(j,:)
             s = s + nstar/(4*PI*norm(node_coordinates(l,:)-node_coordinates(j,:))**3)
             s3 = s3 + intphi_over(l)*cdexp((0,1)*k_wave*norm(node_coordinates(l,:)-node_coordinates(j,:)))/(4*PI*norm(node_coordinates(l,:)-node_coordinates(j,:)))
          end if
       end do
       gauss(j,1) = s
       gauss(j,3) = s3
    end do
  end subroutine setgauss

  double precision function calcsing()
    use omp_lib
    use dbsym
    integer i, nt,iz,ik

    call OMP_SET_NUM_THREADS(4)

    ptr_singular => fsingular3
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call singrate(node_coordinates(i,:),node_coordinates(i,:),i,k_wave,nt)
       gauss(i,4) = sum(jacobian(nt + 1,:,:))
    end do
    !$OMP END PARALLEL DO
  end function calcsing

  subroutine calcomp()
    use omp_lib
    use dbsym
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    ptr_jacobian => fjacobian
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
       intphi_over(i) = folding(i,f,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp

  subroutine calcomp2()
    use omp_lib
    use dbsym
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    ptr_jacobian => fjacobian2
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
       intphi_under(i) = folding(i,f,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp2

  subroutine calcomp3()
    use omp_lib
    use dbsym
    integer i, nt, k1, j

    call OMP_SET_NUM_THREADS(4)

    centres(:) = quadphi_over(:,1)
    weights(:) = quadphi_over(:,2)
    ptr_jacobian => fjacobian
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
       gauss(i,6) = folding(i,f2,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
    do j=2,max_neighbors
       k1 = node_neighbors1(j_tmp,j)
       if (k1 .eq. 0) exit
       call integrate(nstroke_coordinates(k1,:), node_coordinates(k1,:),k1,1)
       gauss(k1, 6) = folding(k1,f3,dim_quad,1)
    end do

  end subroutine calcomp3

  subroutine calcomp4()
    use omp_lib
    use dbsym
    integer i, nt, k1, j

    hval2 = hval
    centres(:) = quadphi_under(:,1)
    weights(:) = quadphi_under(:,2)
    ptr_jacobian => fjacobian2
    call integrate(nstroke_coordinates(j_tmp,:), node_coordinates(j_tmp,:),j_tmp,1)
    do j=2,max_neighbors
       k1 = node_neighbors1(j_tmp,j)
       if (k1 .eq. 0) exit
       i_tmp = k1
       gauss(k1, 6) = gauss(k1, 6) + folding(j_tmp,f4,dim_quad,1)
    end do
    hval2 = hval * hval
    gauss(j_tmp,6) = intphi_under(j_tmp)
    gauss(j_tmp,7) = gauss(j_tmp,4) - (sum(gauss(1:j_tmp-1,6)) + sum(gauss(j_tmp+1:numnodes,6)))/(4*PI)

  end subroutine calcomp4

!       _       _
!      (_)_ __ | |_ ___  __ _
!      | | '_ \| __/ _ \/ _` |
!      | | | | | ||  __/ (_| |
!      |_|_| |_|\__\___|\__, |
!                       |___/

  subroutine singrate(n,z,ip,k_wave,nt)
    use dbsym
    integer, intent(in) :: ip,nt
    double precision, intent(in) :: n(:)
    double precision, intent(in) :: z(:)
    double complex, intent(in) :: k_wave

    double precision rh,ph

    double precision, dimension(dim_3d) :: x, y, p
    double precision, dimension(dim_3d,dim_3d) :: bt

    integer i,nthread

    integer Nk, iz, ik
    double complex jac
    double precision q, ispole
    nthread = nt + 1
    Nk = 4*dim_quad

    p(1) = z(1)
    p(2) = z(2)
    p(3) = z(3)

    ! DONE ugly, 0.5 is value for ellipsoid (*,*,0.5) only
    ispole = 1.0
    if ((p(1)**2 + p(2)**2 .le. 1.0E-12) .and. ((p(3) - axes(3)) .le. 1.0E-12)) then
       p(3) = -axes(3)
       ispole = 0.0
    end if

    do iz=1,dim_quad
       rh = centres(iz)

       do ik=1,Nk
          ph = (2.D0*PI/Nk)*ik
          ! TODO formula (x,y,z)(rh,ph)
          ! do i=1,3
          !    x(i) = (x)s
          !    nodes(nthread,iz,ik,i) = x(i)
          ! end do

          jac =  ptr_singular(axes,p,rh,ph,ispole,k_wave)
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*weights(iz) * jac
       end do
    end do

    return
  end subroutine singrate

  subroutine integrate(n,z,ip,nt)
    use dbsym
    integer, intent(in) :: ip,nt
    double precision, intent(in), dimension(:) :: n
    double precision, intent(in), dimension(:) :: z

    double precision rh,ph

    double precision, dimension(dim_3d) :: x, y
    double precision, dimension(dim_3d,dim_3d) :: bt

    integer m,l,k_ind,i,nthread

    integer Nk, iz, ik
    double precision jac
    nthread = nt + 1
    Nk = 4*dim_quad

    k_ind = MINLOC(n,dim=1)

    do m=1,dim_3d
       do l=1,dim_3d
          bt(m,l) = beta(m,k_ind,l,n,axes)
       end do
    end do

    do iz=1,dim_quad
       rh = centres(iz)

       do ik=1,Nk
          ph = (2.D0*PI/Nk)*ik
          do i=1,3
             x(i) = fx(i,bt,axes,hval2,rh,ph,z)
             nodes(nthread,iz,ik,i) = x(i)
          end do

          jac = dsqrt(ptr_jacobian(bt,axes,hval2,rh,ph,z))
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*weights(iz) * jac

       end do
    end do
    return
  end subroutine integrate

  double complex function folding(ip,f,dim_quad,nt)
    use dbsym
    integer, intent(in) :: ip,nt,dim_quad
    interface
       function f(x,i)
         integer, intent(in) :: i
         double precision, intent(in), dimension(:) :: x
         double precision :: f
       end function f
    end interface

    integer Nk, iz, ik,nthread
    double complex gtmp, tmp

    Nk = 4*dim_quad
    nthread = nt + 1

    gtmp = 0
    do iz=1,dim_quad
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

    Nk = 4*dim_quad
    nthread = nt + 1

    gtmp = (0,0)
    do iz=1,dim_quad
       do ik=1,Nk
          tmp = farr(nthread,iz,ik)
          gtmp = gtmp + real(jacobian(nthread,iz,ik))*tmp
       end do
    end do
    foldingarr = gtmp

  end function foldingarr

end module modinteg
