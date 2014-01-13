module params
  integer :: dim_3d
  integer :: dim_quad
  integer :: max_neighbors
  integer :: numnodes
  double precision :: PI
  double precision :: hval
  double precision :: hval2

  double precision, pointer :: area(:)

  double precision, pointer :: axes(:)

  double precision, pointer :: node_coordinates(:,:)
  double precision, pointer :: normal_coordinates(:,:)
  double precision, pointer :: nstroke_coordinates(:,:)
  double precision, pointer :: intphi_over(:)
  double precision, pointer :: intphi_under(:)

  integer, pointer :: node_neighbors2(:,:)
  integer, pointer :: node_neighbors1(:,:)

  double precision, pointer :: centres(:)
  double precision, pointer :: weights(:)
  double complex, pointer   :: jacobian(:,:,:)
  double precision, pointer :: nodes(:,:,:,:)
  double precision, pointer :: nodes_rho(:,:,:,:)

  double complex :: k_wave
  double precision, pointer   :: sigma(:)
  double complex, pointer     :: q_density(:)

  double complex, pointer     :: gauss(:,:)

  abstract interface
     function iface_f (x,i)
       double precision :: iface_f
       integer, intent(in) :: i
       double precision, intent(in), dimension(:) :: x
     end function iface_f
  end interface
  procedure (iface_f), pointer :: ptr_f => null ()
end module params
