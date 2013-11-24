subroutine set_area(val)
  double precision, intent(in), target :: val(:)
  area => val
end subroutine set_area

subroutine set_nd(val)
  integer, intent(in) :: val
  nd = val
end subroutine set_nd

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

subroutine set_h(val)
  double precision, intent(in) :: val
  h = val
end subroutine set_h

subroutine set_h2(h)
  double precision, intent(in) :: h
  h2 = h
end subroutine set_h2

subroutine set_nz(n)
  integer, intent(in) :: n
  Nz = n
end subroutine set_nz

subroutine set_node_coordinates(points)
  double precision, intent(in), target :: points(:,:)
  node_coordinates => points
end subroutine set_node_coordinates

subroutine set_normal_coordinates(nvectors)
  double precision, intent(in), target :: nvectors(:,:)
  normal_coordinates => nvectors
end subroutine set_normal_coordinates

subroutine set_nstroke_coordinates(nstroke)
  double precision, intent(in), target :: nstroke(:,:)
  nstroke_coordinates => nstroke
end subroutine set_nstroke_coordinates

subroutine set_node_neighbors1(nneigh1)
  integer, intent(in), target :: nneigh1(:,:)
  node_neighbors1 => nneigh1
end subroutine set_node_neighbors1

subroutine set_node_neighbors2(nneigh2)
  integer, intent(in), target :: nneigh2(:,:)
  node_neighbors2 => nneigh2
end subroutine set_node_neighbors2

subroutine set_intphi_over(values)
  double precision, intent(in), target :: values(:)
  intphi_over => values
end subroutine set_intphi_over

subroutine set_intphi_under(values)
  double precision, intent(in), target :: values(:)
  intphi_under => values
end subroutine set_intphi_under

subroutine set_axes(axesi)
  double precision, intent(in), target :: axesi(:)
  axes => axesi
end subroutine set_axes

subroutine set_centres(centresi)
  double precision, intent(in), target :: centresi(:)
  centres => centresi
end subroutine set_centres

subroutine set_c(Ci)
  double precision, intent(in), target :: Ci(:)
  C => Ci
end subroutine set_c

subroutine set_jacobian(jacobiani)
  double complex, intent(in), target :: jacobiani(:,:,:)
  jacobian => jacobiani
end subroutine set_jacobian

subroutine set_nodes(nodesi)
  double precision, intent(in), target :: nodesi(:,:,:,:)
  nodes => nodesi
end subroutine set_nodes

subroutine set_sigma(values)
  double precision, intent(in), target :: values(:)
  sigma => values
end subroutine set_sigma

subroutine set_q(values)
  double complex, intent(in), target :: values(:)
  q => values
end subroutine set_q

subroutine set_k(val)
  double complex, intent(in) :: val
  k = val
end subroutine set_k

subroutine set_gauss(values)
  double complex, intent(in), target :: values(:,:)
  gauss => values
end subroutine set_gauss
