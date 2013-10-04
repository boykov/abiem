module params
  integer :: nd
  integer :: Nz
  integer :: max_neighbors
  integer :: numnodes
  double precision :: PI
  double precision :: h
  double precision :: h2

  double precision, pointer :: area(:)

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

end module params
