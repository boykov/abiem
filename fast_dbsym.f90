module fast_dbsym
  implicit double precision (a-h, o-z)
  double precision, parameter :: PI = 3.14159265358979324D0
contains

  pure double precision function asqrt(x)
    double precision, intent(in) :: x
    asqrt = dsqrt(abs(x))
  end function asqrt

  pure integer function dn(n,m)
    integer, intent(in) :: n,m
    dn = 0
    if (n .eq. m) dn = 1
  end function dn

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

end module fast_dbsym
