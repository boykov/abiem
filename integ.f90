module modinteg
  use params
  use modphi, only : phi
contains

  include 'set_params.f90'
  include 'kernels.f90'

  double precision function f2(x,i)
    use dbsym
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f2 = phi(x,i,hval**2)/norm(x(:) + node_coordinates(i,:) - node_coordinates(1,:))
  end function f2

  double precision function f(x,i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f = phi(x,i,hval**2)
  end function f

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

    ptr_singular => fsingular3
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call singrate(node_coordinates(i,:),node_coordinates(i,:),i,k_wave,centres,weights,nt)
       gauss(i,4) = sum(jacobian(nt + 1,:,:))
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

    ptr_jacobian => fjacobian
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,centres,weights,nt)
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
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,centres,weights,nt)
       intphi_under(i) = folding(i,f,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
  end subroutine calcomp2

  subroutine calcomp3()
    use omp_lib
    use dbsym
    integer i, nt

    call OMP_SET_NUM_THREADS(4)

    ptr_jacobian => fjacobian
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) PRIVATE(nt)
    do i=1,numnodes
       nt = OMP_GET_THREAD_NUM()
       call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,centres,weights,nt)
       gauss(i,6) = folding(i,f2,dim_quad,nt)
    end do
    !$OMP END PARALLEL DO
    gauss(1,6) = 0.0
  end subroutine calcomp3

!       _       _
!      (_)_ __ | |_ ___  __ _
!      | | '_ \| __/ _ \/ _` |
!      | | | | | ||  __/ (_| |
!      |_|_| |_|\__\___|\__, |
!                       |___/

  subroutine singrate(n,z,ip,k_wave,centres,weights,nt)
    use dbsym
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

          jac =  ptr_singular(axes,p,rh,ph,ispole,k_wave)
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*weights(iz) * jac
       end do
    end do

    return
  end subroutine singrate

  subroutine integrate(n,z,ip,centres,weights,nt)
    use dbsym
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

          jac = dsqrt(ptr_jacobian(bt,axes,hval2,rh,ph,z))
          jacobian(nthread,iz,ik) = 2*(PI/Nk)*weights(iz) * jac

       end do
    end do
    return
  end subroutine integrate

  double precision function folding(ip,f,dim_quad,nt)
    use dbsym
    integer, intent(in) :: ip,nt, dim_quad
    interface
       function f(x,i)
         integer, intent(in) :: i
         double precision, intent(in), dimension(:) :: x
         double precision :: f
       end function f
    end interface

    integer Nk, iz, ik,nthread
    double precision gtmp, tmp

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
