module dbsym
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

  pure double precision function fjacobian(bt,axes,h2,rh,ph,z)
    double precision, intent(in) :: bt(:,:)
    double precision, intent(in) :: z(:), axes(:)
    double precision, intent(in) :: h2, rh, ph
    fjacobian = &
         include 'jacobian.f90'
  end function fjacobian

end module dbsym
