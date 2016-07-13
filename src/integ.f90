module modinteg
  use params
  use modphi, only : phi, varphi, deltah
  double precision :: y_tmp(3)
contains

  include 'set_params.f90'
  include 'kernels.f90'

  double complex function fAmn(x,i,j)
    use dbsym
    integer, intent(in) :: i,j
    double precision, intent(in), dimension(:) :: x
    double precision, dimension(size(x)) :: y
    y(:) = x(:) + node_coordinates(i,:) - y_tmp(:)

    fAmn = phi(x,i,hval**2)*Amn(y/2,-y/2,k_wave)
  end function fAmn

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

    call OMP_SET_NUM_THREADS(omp_threads)

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

  double complex function foldingG(N,y,j,k)
    use omp_lib
    use fast_dbsym
    double complex, intent(in) :: k
    double precision, intent(in), dimension(:) :: y
    integer, intent(in) :: j,N
    integer :: l,m, nt
    double precision, dimension(3) :: x
    double complex, dimension(N+1) :: tmp
    double complex, dimension(omp_threads) :: s

    x(:) = y(:) - node_coordinates(j,:)
    tmp = aspherical_hankel_n(N,realpart(k)*norm(x))

    call OMP_SET_NUM_THREADS(omp_threads)

    s(:) = 0.0
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do l = 0,N
       nt = OMP_GET_THREAD_NUM() + 1
       do m = -l, l
          s(nt) = s(nt) + tmp(l+1) * G_harmonic_y(x(:),k,l,m) * intG_x(j,l+1,dim_intG+1+m)
       end do
    end do
    !$OMP END PARALLEL DO
    foldingG = sum(s(:))
  end function foldingG

  subroutine calcomp()
    use omp_lib
    use dbsym
    use fast_dbsym, PI_new => PI, norm_new => norm
    integer :: j_tmp, j_init
    integer i, nt
    integer l,m
    double precision, dimension(3) :: x
    double complex :: tmp

    j_tmp = numnodes
    j_init = numnodes

    if (matrixa6_p) j_init=1

    call OMP_SET_NUM_THREADS(omp_threads)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt, k1, hval2, x)
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

       if (qbx) then
          do iz=1,dim_quad
             do ik=1,4*dim_quad
                x(:) = nodes(nt+1,iz,ik,:)
                cache_phi(nt+1,iz,ik) = phi(x,i,hval**2)
                cache_bessel_jn(nt+1,iz,ik,:) = aspherical_bessel_jn(dim_intG,realpart(k_wave)*norm(x))
             end do
          end do

          do l=0,dim_intG
             do m =-l,l

                do iz=1,dim_quad
                   do ik=1,4*dim_quad
                      x(:) = nodes(nt+1,iz,ik,:)
                      farr(nt+1,iz,ik) = cache_bessel_jn(nt+1,iz,ik,l+1)*G_harmonic_x(x,k_wave,l,m)*cache_phi(nt+1,iz,ik)
                   end do
                end do

                intG_x(i,l+1,dim_intG+1+m) = foldingarr(farr,dim_quad, nt)
             end do
          end do
       end if

       if (gauss6) then
          do j_tmp=j_init,numnodes
             do j=1,max_neighbors
                k1 = node_neighbors2(j_tmp,j)
                if (k1 .eq. i) then
                   int_neighbors2(j_tmp,j) = int_neighbors2(j_tmp,j) + folding(k1,j_tmp,f3,dim_quad,nt)
                   exit
                end if
                if ((k1 .eq. 0)) then
                   exit
                end if
             end do
          end do

          hval2 = hval**2
          call integrate(                &
               fjacobian2,               &
               nstroke_coordinates(i,:), &
               node_coordinates(i,:),    &
               i,                        &
               quadphi_under(:,1),       &
               quadphi_under(:,2),       &
               nt)
          intphi_under(i) = folding(i,i,f,dim_quad,nt)
          do j=1,max_neighbors
             k1 = node_neighbors2(i,j)
             if (k1 .eq. 0) exit
             int_neighbors2(i,j) = int_neighbors2(i,j) + folding(i,k1,f4,dim_quad,nt)
          end do
          hval2 = hval * hval
       end if
    end do
    !$OMP END PARALLEL DO

    if (gauss6) then
       do j_tmp=j_init,numnodes
          gauss(:,6) = 0
          do i=1,numnodes
             do j=1,max_neighbors
                k1 = node_neighbors2(j_tmp,j)
                if (k1 .eq. i) then
                   gauss(k1, 6) = int_neighbors2(j_tmp,j)
                   exit
                end if
                if ((k1 .eq. 0)) then
                   if (qbx_gauss6) then
                      gauss(i,6) = gauss(i, 6) + foldingG(dim_intG,node_coordinates(j_tmp,:),i,k_wave)*(4*PI)
                      else
                      call integrate(                &
                           fjacobian,                &
                           nstroke_coordinates(i,:), &
                           node_coordinates(i,:),    &
                           i,                        &
                           quadphi_over(:,1),        &
                           quadphi_over(:,2),        &
                           0)
                      gauss(i,6) = gauss(i, 6) + folding(i,j_tmp,f2,dim_quad,0)
                   end if
                   exit
                end if
             end do
          end do

          gauss(j_tmp,7) = - (sum(gauss(1:j_tmp-1,6)) + sum(gauss(j_tmp+1:numnodes,6)))/(4*PI)
       end do
    end if

  end subroutine calcomp

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

  double complex function foldingarr(farr,dim_quad,nt)
    integer, intent(in) :: nt, dim_quad
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
