module params
  integer :: nd
  integer :: Nz
  integer :: max_neighbors
  integer :: numnodes
  double precision :: PI
  double precision :: h
  double precision :: h2

  double precision :: area

  double precision, pointer :: axes(:)

  double precision, pointer :: node_coordinates(:,:)
  double precision, pointer :: normal_coordinates(:,:)
  double precision, pointer :: nstroke_coordinates(:,:)
  double precision, pointer :: intphi_over(:)

  integer, pointer :: node_neighbors2(:,:)
  integer, pointer :: node_neighbors1(:,:)

  double precision, pointer :: centres(:)
  double precision, pointer :: C(:)
  double complex, pointer   :: jacobian(:,:,:)
  double precision, pointer :: nodes(:,:,:,:)

contains

  double precision function asqrt(x)
    double precision, intent(in) :: x
    asqrt = dsqrt(abs(x))
  end function asqrt

  integer function dn(n,m)
    integer :: n,m
    dn = 0
    if (n .eq. m) dn = 1
  end function dn
end module params
