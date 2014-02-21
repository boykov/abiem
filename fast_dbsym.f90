module fast_dbsym
  implicit double precision (a-h, o-z)
  double precision, parameter :: PI = 3.14159265358979324D0
contains
  include 'msphj.for'

  pure double precision function asqrt(x)
    double precision, intent(in) :: x
    asqrt = dsqrt(abs(x))
  end function asqrt

  pure integer function dn(n,m)
    integer, intent(in) :: n,m
    dn = 0
    if (n .eq. m) dn = 1
  end function dn

  double complex function cderf(z)
    double complex, intent(in) :: z
    cderf = &
         include 'cderf.f90'
  end function cderf

  double complex function A_s(x,y,sigm,k)
    double precision, intent(in) :: x(:),y(:)
    double precision, intent(in) :: sigm
    double complex, intent(in) :: k
    double precision :: v,rh
    rh = norm(x-y)
    A_s = &
         include 'A_s.f90'
  end function A_s

  real(16) function spherical_bessel_jn(n, x)
    integer, intent(in) :: n
    real(16), intent(in) :: x
    integer :: nm
    real(16) :: sj(0:n), dj(0:n)
    call sphj(n, x, nm, sj, dj)
    if (nm /= n) print *, "spherical_bessel_jn: sphj didn't converge"
    spherical_bessel_jn = sj(n)
  end function spherical_bessel_jn

  real(16) function spherical_bessel_yn(n, x)
    integer, intent(in) :: n
    real(16), intent(in) :: x
    integer :: nm
    real(16) :: sj(0:n), dj(0:n)
    call sphy(n, x, nm, sj, dj)
    if (nm /= n) print *, "spherical_bessel_jn: sphj didn't converge"
    spherical_bessel_yn = sj(n)
  end function spherical_bessel_yn

  double complex function spherical_hankel_n(n,x)
    double precision, intent(in) :: x
    integer, intent(in) :: n
    spherical_hankel_n = spherical_bessel_jn(n,x) + (0,1)*spherical_bessel_yn(n,x)
  end function spherical_hankel_n

  double complex function spherical_hankel(n,x)
    double precision, intent(in) :: x
    integer, intent(in) :: n
    spherical_hankel = spherical_bessel_j(n,x) + (0,1)*spherical_bessel_y(n,x)
  end function spherical_hankel

  double precision function spherical_bessel_y(n,x)
    double precision, intent(in) :: x
    integer, intent(in) :: n
    spherical_bessel_y = &
         include 'spherical_bessel_y.f90'
  end function spherical_bessel_y

  double precision function spherical_bessel_j(n,x)
    double precision, intent(in) :: x
    integer, intent(in) :: n
    spherical_bessel_j = &
         include 'spherical_bessel_j.f90'
  end function spherical_bessel_j

  double complex function spherical_harmonic(l,m,theta,phi)
    double precision, intent(in) :: theta,phi
    integer, intent(in) :: l,m
    spherical_harmonic = &
         include 'spherical_harmonic.f90'
  end function spherical_harmonic

  double complex function spherical_harmonic_(l,m,theta,phi)
    double precision, intent(in) :: theta,phi
    integer, intent(in) :: l,m

    spherical_harmonic_ = spherical_harmonic(l,m,phi,theta)/cdexp((0,1)*m*theta)
  end function spherical_harmonic_

  double complex function G_y(y,k,l,m)
    double precision, intent(in) :: y(:)
    double complex, intent(in) :: k
    integer, intent(in) :: l,m
    double precision :: ry, phiy, thetay

    ry = norm(y)
    phiy = dacos(y(3)/ry)
    thetay = datan(y(2)/y(1))

    G_y = (0,1)*realpart(k) * spherical_hankel_n(l, realpart(k)*ry) * &
         cdexp((0,1)*(-thetay)*m)*spherical_harmonic_(l,m,thetay,phiy)

  end function G_y

  double complex function G_x(x,k,l,m)
    double precision, intent(in) :: x(:)
    double complex, intent(in) :: k
    integer, intent(in) :: l,m
    double precision :: rx, phix, thetax

    rx = norm(x)
    phix = dacos(x(3)/rx)
    thetax = datan(x(2)/x(1))

    G_x =  spherical_bessel_jn(l, realpart(k)*rx) * &
         cdexp((0,1)*thetax*m)*spherical_harmonic_(l,m,thetax,phix)

  end function G_x

  ! formula in [[file:~/downloads/pub/papers/epstein2012convergence.pdf][epstein]]

  double complex function G_2(x,y,k)
    double precision, intent(in) :: x(:),y(:)
    double complex, intent(in) :: k
    integer :: l,m
    double complex :: s

    s = 0.0
    do l = 0,20
       do m = -l, l
          s = s + G_y(y,k,l,m)* G_x(x,k,l,m)
       end do
    end do

    G_2 = s
  end function G_2

  double complex function G_(x,y,k)
    double precision, intent(in) :: x(:),y(:)
    double complex, intent(in) :: k
    double precision :: rx, phix, thetax
    double precision :: ry, phiy, thetay
    double complex :: s, s2
    integer :: l,m

    rx = norm(x)
    ry = norm(y)
    phix = dacos(x(3)/rx)
    phiy = dacos(y(3)/ry)
    thetax = datan(x(2)/x(1))
    thetay = datan(y(2)/y(1))

    s = 0.0
    do l = 0,20
       s2 = 0.0
       do m = -l, l
          s2 = s2 + cdexp((0,1)*(thetax-thetay)*m) * &
               spherical_harmonic_(l,m,thetax,phix)* &
               spherical_harmonic_(l,m,thetay,phiy)
       end do
       s = s + (0,1)*realpart(k)* spherical_bessel_jn(l, realpart(k)*rx) * &
            spherical_hankel_n(l, realpart(k)*ry) * s2
    end do
    G_ = s
  end function G_

  double precision function foo2()
    foo2 = 5.0
  end function foo2

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
    if (dabs(dble(y)) .ge. 3) then
       alpha = dabs(dble(y/1e+5))
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
       wmod = dexp(dble(-y**2))*w(z)
    end if
  end function wmod

  double complex function w(z)
    use toms
    real(dp) x,y,u,v
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

end module fast_dbsym
