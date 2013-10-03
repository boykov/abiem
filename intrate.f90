module modint
    use params
contains
  subroutine integrate(n,z,ip,nt)
    
    integer, intent(in) :: ip,nt
    double precision, intent(in), dimension(nd) :: n
    double precision, intent(in), dimension(nd) :: z

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
          bt(m,l) = &
               include 'beta.f90'
       end do
    end do

    do iz=1,Nz
       rh = centres(iz)

       do ik=1,Nk
          ph = (2.D0*PI/Nk)*ik
          do i=1,3
             x(i) = &
                  include 'x.f90'
             nodes(nthread,iz,ik,i) = x(i)
          end do

          jac = dsqrt (&
               include 'jacobian.f90'
          )
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
    double precision, dimension(nd) :: getnode
    getnode(:) = nodes(nt,iz,ik,:)
  end function getnode

  double complex function foldingarr(ip,farr,nt)
    integer, intent(in) :: ip,nt
    double complex, dimension(4,Nz,4*Nz) :: farr
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
end module modint

