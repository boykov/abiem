function ylm (l, m, thrad, phirad)
  implicit none
  !--------------------------------------------------------------------
  ! Computes the spherical harmonic Y_lm (theta,phi) using the
  ! reduced rotation matrix d^l_{m 0} (theta) and using the
  ! external function fac10(n) = factorial(n)/10**n
  !--------------------------------------------------------------------
  ! input: angular momentum quantum numbers l, m (integers)
  !        angles theta and phi (radian)
  ! -------------------------------------------------------------------
  ! Reference: D.M. Brink and G.R. Satchler, Angular Momentum,
  !            second edition, Oxford University Press, p.22 and p. 145
  ! -------------------------------------------------------------------
  integer, parameter  :: wp = kind(1.0d0)  ! working precision = double (portable)
  !--------------------------------------------------------------------
  !   local constants
  !--------------------------------------------------------------------
  real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
  complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
  !--------------------------------------------------------------------
  !   formal arguments
  !--------------------------------------------------------------------
  integer, intent(in)  :: l, m
  real(wp), intent(in) :: thrad, phirad
  !--------------------------------------------------------------------
  !   local variables
  !--------------------------------------------------------------------
  integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, it, iphase, &
       ia, ib, ic
  real(wp) :: sqrt_fac, sumt, denom, term, dlm0, const, cosb2, sinb2
  complex(wp) :: ylm, exphi
  !--------------------------------------------------------------------
  !   external function
  !--------------------------------------------------------------------
  ! real(wp), external :: fac10
  !--------------------------------------------------------------------
  !  program starts here
  !  first calculate d^l_{m 0} (theta)
  !--------------------------------------------------------------------
  cosb2 = cos(thrad/2.0_wp)
  sinb2 = sin(thrad/2.0_wp)
  !--------------------------------------------------------------------
  ! determine lower and upper limits for summation index it; these
  ! are derived from the requirement that all factorials n! in the
  ! denominator are restricted to values with n >=0.
  !--------------------------------------------------------------------
  itmin1 = 0
  itmin2 = m
  itmin = max(itmin1,itmin2)
  itmax1 = l+m
  itmax2 = l
  itmax = min(itmax1,itmax2)
  !  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
  sqrt_fac = sqrt( fac10(l+m) * fac10(l-m) * fac10(l) * fac10(l) )
  !
  sumt = 0.0_wp
  do it = itmin, itmax
     iphase = (-1)**it
     ia = l + m - it
     ib = l - it
     ic = it - m
     !     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
     denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
     term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
     sumt = sumt + term
  end do
  dlm0 = sqrt_fac * sumt
  !--------------------------------------------------------------------
  !  now compute Y_{l m} (theta,phi) from d^l_{m 0} (theta)
  !--------------------------------------------------------------------
  const = sqrt( (2.0_wp *l + 1.0_wp) / (4.0_wp * pi) )
  exphi = exp( eye * m * phirad )
  ylm = const * exphi * dlm0
  !
  return
end function ylm


function fac10 (n)
  implicit none
  ! -----------------------------------------------
  ! function fac10(n) calculates factorial(n)/10**n
  ! -----------------------------------------------
  ! input: integer n >= 0 (you may want to check this
  !        in the program calling this function)
  ! -----------------------------------------------
  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
  !------------------------------------------------
  !      formal arguments
  !------------------------------------------------
  integer, intent(in) :: n
  !------------------------------------------------
  !      local variables
  !------------------------------------------------
  integer :: i
  real(wp) :: fac10, q
  ! -----------------------------------------------
  if (n == 0) then
     fac10 = 1.0_wp
  else
     fac10 = 1.0_wp
     q = 1.0_wp
     do i = 1, n
        fac10 = fac10 * q / 10.0_wp
        q = q + 1.0_wp
     end do
  endif
  !
  return
end function fac10
