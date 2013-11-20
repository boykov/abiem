module fast_dbsym
  implicit double precision (a-h, o-z)
  double precision, parameter :: PI = 3.14159265358979324D0
contains

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
