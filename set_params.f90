subroutine set_dp2d_ptr(name, val)
  character(128) :: name
  double precision, intent(in), target :: val(:,:)
  if (name .eq. "node_coordinates") node_coordinates => val
  if (name .eq. "normal_coordinates") normal_coordinates => val
  if (name .eq. "nstroke_coordinates") nstroke_coordinates => val
  if (name .eq. "quadphi_over") quadphi_over => val
  if (name .eq. "quadphi_under") quadphi_under => val
end subroutine set_dp2d_ptr

subroutine set_dp1d_ptr(name, val)
  character(128) :: name
  double precision, intent(in), target :: val(:)
  if (name .eq. "centres") centres => val
  if (name .eq. "weights") weights => val
  if (name .eq. "area") area => val
  if (name .eq. "axes") axes => val
  if (name .eq. "intphi_over") intphi_over => val
  if (name .eq. "intphi_under") intphi_under => val
  if (name .eq. "sigma") sigma => val
end subroutine set_dp1d_ptr

subroutine set_dp_ptr(name, val)
  character(128) :: name
  double precision, intent(in) :: val
  if (name .eq. "PI") PI = val
  if (name .eq. "hval") hval = val
  if (name .eq. "hval2") hval2 = val
end subroutine set_dp_ptr

subroutine set_i_ptr(name, val)
  character(128) :: name
  integer, intent(in) :: val
  if (name .eq. "dim_3d") dim_3d = val
  if (name .eq. "dim_quad") dim_quad = val
  if (name .eq. "max_neighbors") max_neighbors = val
  if (name .eq. "numnodes") numnodes = val
end subroutine set_i_ptr

subroutine set_area(val)
  double precision, intent(in), target :: val(:)
  area => val
end subroutine set_area

subroutine set_dim_3d(val)
  integer, intent(in) :: val
  dim_3d = val
end subroutine set_dim_3d

subroutine set_pi(val)
  double precision, intent(in) :: val
  PI = val
end subroutine set_pi

subroutine set_max_neighbors(val)
  integer, intent(in) :: val
  max_neighbors = val
end subroutine set_max_neighbors

subroutine set_numnodes(val)
  integer, intent(in) :: val
  numnodes = val
end subroutine set_numnodes

subroutine set_hval(val)
  double precision, intent(in) :: val
  hval = val
end subroutine set_hval

subroutine set_hval2(val)
  double precision, intent(in) :: val
  hval2 = val
end subroutine set_hval2

subroutine set_dim_quad(val)
  integer, intent(in) :: val
  dim_quad = val
end subroutine set_dim_quad

subroutine set_node_coordinates(val)
  double precision, intent(in), target :: val(:,:)
  node_coordinates => val
end subroutine set_node_coordinates

subroutine set_normal_coordinates(val)
  double precision, intent(in), target :: val(:,:)
  normal_coordinates => val
end subroutine set_normal_coordinates

subroutine set_nstroke_coordinates(val)
  double precision, intent(in), target :: val(:,:)
  nstroke_coordinates => val
end subroutine set_nstroke_coordinates

subroutine set_node_neighbors1(val)
  integer, intent(in), target :: val(:,:)
  node_neighbors1 => val
end subroutine set_node_neighbors1

subroutine set_node_neighbors2(val)
  integer, intent(in), target :: val(:,:)
  node_neighbors2 => val
end subroutine set_node_neighbors2

subroutine set_intphi_over(val)
  double precision, intent(in), target :: val(:)
  intphi_over => val
end subroutine set_intphi_over

subroutine set_intphi_under(val)
  double precision, intent(in), target :: val(:)
  intphi_under => val
end subroutine set_intphi_under

subroutine set_axes(val)
  double precision, intent(in), target :: val(:)
  axes => val
end subroutine set_axes

subroutine set_centres(val)
  double precision, intent(in), target :: val(:)
  centres => val
end subroutine set_centres

subroutine set_weights(val)
  double precision, intent(in), target :: val(:)
  weights => val
end subroutine set_weights

subroutine set_jacobian(val)
  double complex, intent(in), target :: val(:,:,:)
  jacobian => val
end subroutine set_jacobian

subroutine set_nodes(val)
  double precision, intent(in), target :: val(:,:,:,:)
  nodes => val
end subroutine set_nodes

subroutine set_sigma(val)
  double precision, intent(in), target :: val(:)
  sigma => val
end subroutine set_sigma

subroutine set_q_density(val)
  double complex, intent(in), target :: val(:)
  q_density => val
end subroutine set_q_density

subroutine set_k_wave(val)
  double complex, intent(in) :: val
  k_wave = val
end subroutine set_k_wave

subroutine set_gauss(val)
  double complex, intent(in), target :: val(:,:)
  gauss => val
end subroutine set_gauss
