module intphi
  use params
  use modint
  use modphi, only : phi

contains
  
  subroutine inomp(i,nt)
    integer, intent(in) :: i, nt
    call integrate(nstroke_coordinates(i,:),node_coordinates(i,:),i,nt)
    intphi_over(i) = folding(i,f,nt)
  end subroutine inomp

  double precision function f(x,i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    f = phi(x,i,h**2) + 0
  end function f
  
end module intphi
