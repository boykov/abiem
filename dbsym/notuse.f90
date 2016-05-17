
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

