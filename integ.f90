module modinteg
  use params
  use modphi, only : phi, varphi, deltah
contains

  include 'set_params.f90'
  include 'kernels.f90'

  double complex function f2(x,i,j)
    use dbsym
    integer, intent(in) :: i,j
    double precision, intent(in), dimension(:) :: x
    double precision, dimension(size(x)) :: y
    y(:) = x(:) + node_coordinates(i,:) - node_coordinates(j,:)

    f2 = cdexp((0,1)*k_wave*norm(y))*phi(x,i,hval**2)/norm(y)
  end function f2

  double complex function one(x,i,j)
    integer, intent(in) :: i,j
    double precision, intent(in), dimension(:) :: x
    one = 1
  end function one

  double complex function f(x,i,j)
    integer, intent(in) :: i,j
    double precision, intent(in), dimension(:) :: x
    f = phi(x,i,hval**2)
  end function f

  double complex function f3(x,i,j)
    use dbsym
    integer, intent(in) :: i,j
    double precision, intent(in), dimension(:) :: x
    double precision, dimension(size(x)) :: y
    y(:) = x(:) + node_coordinates(i,:) - node_coordinates(j,:)

    f3 = cdexp((0,1)*k_wave*norm(y))*(1 - varphi(y,hval**2))*phi(x,i,hval**2)/norm(y)
  end function f3

  double complex function f4(x,i,j)
    use dbsym
    integer, intent(in) :: i,j
    double precision, intent(in), dimension(:) :: x
    double precision, dimension(size(x)) :: y
    y(:) = x(:) + node_coordinates(i,:) - node_coordinates(j,:)

    f4 = cdexp((0,1)*k_wave*norm(x))*(varphi(x,hval**2))*phi(y,j,hval**2)
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

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(s,s3,x,y)
    do j=1,numnodes
       s = 0.0
       s2 = 0.0
       s3 = 0.0
       do l=1,numnodes
          if (l .ne. j) then
             sigm = sigmaij(l,j)
             x = node_coordinates(l,:)
             y = node_coordinates(j,:)
             nstar = sum(normal_coordinates(l,:)*(x-y)*intphi_over(l))
             s = s + nstar*(-Bmn(x,y,DCMPLX(0,0)))
             s3 = s3 + intphi_over(l)*Amn(x,y,k_wave)
          end if
       end do
       gauss(j,1) = s
       gauss(j,3) = s3
    end do
    !$OMP END PARALLEL DO
  end subroutine setgauss

  double precision function calcsing()
    use omp_lib
    use dbsym
    integer i, nt,iz,ik

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call singrate(              &
            fsingular3,            &
            node_coordinates(i,:), &
            node_coordinates(i,:), &
            i,                     &
            k_wave,                &
            quadsingular(:,1),     &
            quadsingular(:,2),     &
            nt)
       gauss(i,4) = folding(i,i,one,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
  end function calcsing

  subroutine setup_calcomp ()
    use dbsym
    ptr_jacobian => fjacobian
  end subroutine setup_calcomp

  subroutine calcomp()
    use omp_lib
    use dbsym
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(                &
            fjacobian,                &
            nstroke_coordinates(i,:), &
            node_coordinates(i,:),    &
            i,                        &
            quadphi_over(:,1),        &
            quadphi_over(:,2),        &
            nt)
       intphi_over(i) = folding(i,i,f,dim_quad,nt)

       call integrate(                &
            fjacobian2,               &
            nstroke_coordinates(i,:), &
            node_coordinates(i,:),    &
            i,                        &
            quadphi_under(:,1),       &
            quadphi_under(:,2),       &
            nt)
       intphi_under(i) = folding(i,i,f,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp

  subroutine calcomp2()
    use omp_lib
    use dbsym
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(                &
            fjacobian2,               &
            nstroke_coordinates(i,:), &
            node_coordinates(i,:),    &
            i,                        &
            quadphi_under(:,1),       &
            quadphi_under(:,2),       &
            nt)
       intphi_under(i) = folding(i,i,f,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp2

  subroutine calcomp3(j_tmp)
    use omp_lib
    use dbsym
    integer, intent(in) :: j_tmp
    integer i, nt, k1, j

    call OMP_SET_NUM_THREADS(4)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(                &
            fjacobian,                &
            nstroke_coordinates(i,:), &
            node_coordinates(i,:),    &
            i,                        &
            quadphi_over(:,1),        &
            quadphi_over(:,2),        &
            nt)
       gauss(i,6) = folding(i,j_tmp,f2,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
    do j=1,max_neighbors
       k1 = node_neighbors1(j_tmp,j)
       if (k1 .eq. 0) exit
       call integrate(                 &
            fjacobian,                 &
            nstroke_coordinates(k1,:), &
            node_coordinates(k1,:),    &
            k1,                        &
            quadphi_over(:,1),         &
            quadphi_over(:,2),         &
            1)
       gauss(k1, 6) = folding(k1,j_tmp,f3,dim_quad,1)
    end do

  end subroutine calcomp3

  subroutine calcomp4(j_tmp)
    use omp_lib
    use dbsym
    integer, intent(in) :: j_tmp
    integer i, nt, k1, j

    hval2 = hval**2
    call integrate(                    &
         fjacobian2,                   &
         nstroke_coordinates(j_tmp,:), &
         node_coordinates(j_tmp,:),    &
         j_tmp,                        &
         quadphi_under(:,1),           &
         quadphi_under(:,2),           &
         1)
    do j=1,max_neighbors
       k1 = node_neighbors1(j_tmp,j)
       if (k1 .eq. 0) exit
       gauss(k1, 6) = gauss(k1, 6) + folding(j_tmp,k1,f4,dim_quad,1)
    end do
    hval2 = hval * hval

    gauss(j_tmp,7) = gauss(j_tmp,4) - (sum(gauss(1:j_tmp-1,6)) + sum(gauss(j_tmp+1:numnodes,6)))/(4*PI)

  end subroutine calcomp4

!       _       _
!      (_)_ __ | |_ ___  __ _
!      | | '_ \| __/ _ \/ _` |
!      | | | | | ||  __/ (_| |
!      |_|_| |_|\__\___|\__, |
!                       |___/

  subroutine singrate(f,n,z,ip,k_wave,centres,weights,nt)
    use dbsym
    interface
       function f(axes,p,rh,ph,ispole,k)
         double complex :: f
         double precision, intent(in) :: axes(:), p(:)
         double precision, intent(in) :: rh, ph, ispole
         double complex, intent(in) :: k
       end function f
    end interface

    integer, intent(in) :: ip,nt
    double precision, intent(in) :: n(:)
    double precision, intent(in) :: z(:)
    double precision, intent(in) :: centres(:)
    double precision, intent(in) :: weights(:)
    double complex, intent(in) :: k_wave

    double precision rh,ph

    double precision, dimension(dim_3d) :: x, y, p
    double precision, dimension(dim_3d,dim_3d) :: bt

    integer i,nthread

    integer Nk, iz, ik, dim_quad
    double complex jac
    double precision q, ispole
    nthread = nt + 1
    dim_quad = size(centres,1)
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

          jac =  f(axes,p,rh,ph,ispole,k_wave)
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*weights(iz) * jac
       end do
    end do

    return
  end subroutine singrate

  subroutine integrate(f,n,z,ip,centres,weights,nt)
    use dbsym
    interface
       function f (bt,axes,h2,rh,ph,z)
         double precision :: f
         double precision, intent(in) :: bt(:,:)
         double precision, intent(in) :: z(:), axes(:)
         double precision, intent(in) :: h2, rh, ph
       end function f
    end interface

    integer, intent(in) :: ip,nt
    double precision, intent(in), dimension(:) :: n
    double precision, intent(in), dimension(:) :: z
    double precision, intent(in), dimension(:) :: centres
    double precision, intent(in), dimension(:) :: weights

    double precision rh,ph

    double precision, dimension(dim_3d) :: x, y
    double precision, dimension(dim_3d,dim_3d) :: bt

    integer m,l,k_ind,i,nthread

    integer Nk, iz, ik, dim_quad
    double precision jac
    nthread = nt + 1
    dim_quad = size(centres,1)
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

          jac = dsqrt(f(bt,axes,hval2,rh,ph,z))
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*weights(iz) * jac

       end do
    end do
    return
  end subroutine integrate

  double complex function folding(ip,jp,f,dim_quad,nt)
    use dbsym
    integer, intent(in) :: ip,jp,nt, dim_quad
    interface
       function f(x,i,j)
         integer, intent(in) :: i,j
         double precision, intent(in), dimension(:) :: x
         double complex :: f
       end function f
    end interface

    integer Nk, iz, ik,nthread
    double complex gtmp, tmp

    Nk = 4*dim_quad
    nthread = nt + 1

    gtmp = 0
    do iz=1,dim_quad
       do ik=1,Nk
          tmp = f(nodes(nthread,iz,ik,:),ip,jp)
          gtmp = gtmp + (jacobian(nthread,iz,ik))*tmp
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
          gtmp = gtmp + (jacobian(nthread,iz,ik))*tmp
       end do
    end do
    foldingarr = gtmp

  end function foldingarr

end module modinteg
