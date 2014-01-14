subroutine set_dp2d_ptr(name, val)
  character(128) :: name
  double precision, intent(in), target :: val(:,:)
  if (name .eq. "node_coordinates") node_coordinates => val
  if (name .eq. "normal_coordinates") normal_coordinates => val
  if (name .eq. "nstroke_coordinates") nstroke_coordinates => val
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
  if (name .eq. "counter") counter => val
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

subroutine set_i2d_ptr(name, val)
  character(128) :: name
  integer, intent(in), target :: val(:,:)
  if (name .eq. "node_neighbors1") node_neighbors1 => val
  if (name .eq. "node_neighbors2") node_neighbors2 => val
end subroutine set_i2d_ptr

subroutine set_dc3d_ptr(name, val)
  character(128) :: name
  double complex, intent(in), target :: val(:,:,:)
  if (name .eq. "jacobian") jacobian => val
end subroutine set_dc3d_ptr  

subroutine set_dp4d_ptr(name, val)
  character(128) :: name
  double precision, intent(in), target :: val(:,:,:,:)
  if (name .eq. "nodes") nodes => val
end subroutine set_dp4d_ptr

subroutine set_dc1d_ptr(name, val)
  character(128) :: name
  double complex, intent(in), target :: val(:)
  if (name .eq. "q_density") q_density => val
end subroutine set_dc1d_ptr  

subroutine set_dc2d_ptr(name, val)
  character(128) :: name
  double complex, intent(in), target :: val(:,:)
  if (name .eq. "gauss") gauss => val
end subroutine set_dc2d_ptr  

subroutine set_dc_ptr(name, val)
  character(128) :: name
  double complex, intent(in) :: val
  if (name .eq. "k_wave") k_wave = val
end subroutine set_dc_ptr  
