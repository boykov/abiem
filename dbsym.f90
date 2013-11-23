module dbsym
  implicit double precision (a-h, o-z)
  double precision, parameter :: PI = 3.14159265358979324D0
  abstract interface
     function ifjacobian (bt,axes,h2,rh,ph,z)
       double precision :: ifjacobian
       double precision, intent(in) :: bt(:,:)
       double precision, intent(in) :: z(:), axes(:)
       double precision, intent(in) :: h2, rh, ph
     end function ifjacobian
  end interface
  procedure (ifjacobian), pointer :: ptr_jacobian => null ()
  abstract interface
     function ifsingular(axes,p,rh,ph,ispole,k)
       double complex :: ifsingular
       double precision, intent(in) :: axes(:), p(:)
       double precision, intent(in) :: rh, ph, ispole
       double complex, intent(in) :: k
     end function ifsingular
  end interface
  procedure (ifsingular), pointer :: ptr_singular => null ()
contains
  pure double precision function norm2(v)
    double precision, intent(in) :: v(:)
    norm2 = sum(v(:)*v(:))
  end function norm2

  pure double precision function norm(v)

    double precision, intent(in) :: v(:)
    norm = dsqrt(norm2(v))
  end function norm

  double complex function wmod(z)
    use toms
    real(dp) x,y,u,v
    double complex, intent(in) :: z
    real(dp) alpha
    integer(kind=16) m,q,jm
    logical flag

    x = REAL(z)
    y = AIMAG(z)
    if (dabs(y) .ge. 3) then
       alpha = dabs(y/1e+5)
       m = 1e+20
       jm = y/alpha
       q = (real(2*m,dp)/PI)*alpha*x
       call wofz_mod(alpha,m,q,jm,u,v,flag)
       if (flag) then
          ! print *,"Alert wofz_mod"
          ! stop
       end if
       wmod = DCMPLX(u,v)
    else
       wmod = dexp(-y**2)*w(z)
    end if
  end function wmod

  double complex function w(z)
    use toms
    double precision x,y,u,v
    double complex, intent(in) :: z
    logical flag

    x = REAL(z)
    y = AIMAG(z)
    call WOFZ(x,y,u,v,flag)
    if (flag) then
       ! print *,"Alert wofz"
       ! stop
    end if
    w = DCMPLX(u,v)
  end function w

  pure double precision function asqrt(x)
    double precision, intent(in) :: x
    asqrt = dsqrt(abs(x))
  end function asqrt

  pure integer function dn(n,m)
    integer, intent(in) :: n,m
    dn = 0
    if (n .eq. m) dn = 1
  end function dn

  pure double precision function beta(m,k,l,n,axes)
    integer, intent(in) :: m,k,l
    double precision, intent(in) :: n(:), axes(:)
    beta = &
    include 'beta.f90'
  end function beta

  pure double precision function fx(i,bt,axes,h2,rh,ph,z)
    integer, intent(in) :: i
    double precision, intent(in) :: bt(:,:)
    double precision, intent(in) :: z(:), axes(:)
    double precision, intent(in) :: h2, rh, ph
    fx = &
         include 'x.f90'
  end function fx

  pure double complex function fsingular(axes,p,rh,ph,ispole,k)
    double precision, intent(in) :: axes(:), p(:)
    double precision, intent(in) :: rh, ph, ispole
    double complex, intent(in) :: k
    fsingular = &
         include 'singular.f90'
  end function fsingular

  pure double complex function fsingular2(axes,p,rh,ph,ispole,k)
    double precision, intent(in) :: axes(:), p(:)
    double precision, intent(in) :: rh, ph, ispole
    double complex, intent(in) :: k
    fsingular2 = &
         include 'singular2.f90'
  end function fsingular2

  pure double precision function fjacobian(bt,axes,h2,rh,ph,z)
    double precision, intent(in) :: bt(:,:)
    double precision, intent(in) :: z(:), axes(:)
    double precision, intent(in) :: h2, rh, ph
    fjacobian = &
         include 'jacobian.f90'
  end function fjacobian

  pure double precision function fjacobian2(bt,axes,h2,rh,ph,z)
    double precision, intent(in) :: bt(:,:)
    double precision, intent(in) :: z(:), axes(:)
    double precision, intent(in) :: h2, rh, ph
    fjacobian2 = &
         include 'jacobian2.f90'
  end function fjacobian2

  double complex function Amn(x,y,k)
    double precision, intent(in) :: x(:),y(:)
    double complex, intent(in) :: k
    double precision :: v,rh
    rh = norm(x-y)
    Amn = &
         include 'Amn.f90'
  end function Amn

  double complex function Bmn(x,y,k)
    double precision, intent(in) :: x(:),y(:)
    double complex, intent(in) :: k
    double precision :: v,rh
    rh = norm(x-y)
    Bmn = &
         include 'Bmn.f90'
  end function Bmn

  double complex function limdA(sigm,k)
    double precision, intent(in) :: sigm
    double complex, intent(in) :: k
    limdA = &
         include 'limdA.f90'
  end function limdA

  double complex function dA(x,y,sigm,k)
    double precision, intent(in) :: x(:),y(:)
    double precision, intent(in) :: sigm
    double complex, intent(in) :: k
    double precision :: v,rh
    rh = norm(x-y)
    dA = &
         include 'dA.f90'
  end function dA

  double complex function limA(sigm,k)
    double precision, intent(in) :: sigm
    double complex, intent(in) :: k
    limA = &
         include 'limA.f90'
  end function limA

  double complex function A(x,y,sigm,k)
    double precision, intent(in) :: x(:),y(:)
    double precision, intent(in) :: sigm
    double complex, intent(in) :: k
    double precision :: v,rh
    rh = norm(x-y)
    A = &
         include 'A.f90'
  end function A

end module dbsym
