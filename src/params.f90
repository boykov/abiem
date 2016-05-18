module params
  double complex, pointer :: int_neighbors2(:,:)
  double complex, pointer :: int_neighbors1(:,:)
  integer :: numnodes
  double complex, pointer :: valG_y(:,:,:)
  integer :: dim_quad
  logical :: use_int_neighbors_p
  double precision, pointer :: area(:)
  integer :: dim_3d
  double precision, pointer :: quadphi_over(:,:)
  double precision, pointer :: nstroke_coordinates(:,:)
  logical :: qbx_gauss6
  double complex, pointer :: intG_x(:,:,:)
  integer :: max_neighbors
  double precision, pointer :: cache_bessel_jn(:,:,:,:)
  double complex, pointer :: gauss(:,:)
  integer, pointer :: node_neighbors2(:,:)
  integer, pointer :: node_neighbors1(:,:)
  double precision, pointer :: nodes(:,:,:,:)
  double precision, pointer :: intphi_under(:)
  integer :: dim_intG
  double precision, pointer :: node_coordinates(:,:)
  double precision, pointer :: quadphi_under(:,:)
  logical :: matrixa6_p
  double precision :: hval
  logical :: gauss6
  double complex, pointer :: q_density(:)
  double precision, pointer :: axes(:)
  double precision, pointer :: intphi_over(:)
  double complex :: k_wave
  double complex, pointer :: jacobian(:,:,:)
  double precision :: PI
  double precision, pointer :: normal_coordinates(:,:)
  double precision :: hval2
  double complex, pointer :: farr(:,:,:)
  double precision, pointer :: quadsingular(:,:)
  logical :: qbx
  double precision, pointer :: sigma(:)
  double precision, pointer :: counter(:)
  double complex, pointer :: cache_phi(:,:,:)
  double precision, pointer :: weights(:)
  double precision, pointer :: centres(:)
end module params