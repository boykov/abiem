module modphi
  integer :: max_neighbors
  integer :: numnodes
  double precision :: h

  double precision, pointer :: node_coordinates(:,:)
  double precision, pointer :: normal_coordinates(:,:)
  double precision, pointer :: nstroke_coordinates(:,:)

  double precision, pointer :: axes(:)

  integer, pointer :: node_neighbors2(:,:)
  integer, pointer :: node_neighbors1(:,:)


contains

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

  subroutine set_axes(axesi)
    double precision, intent(in), target :: axesi(:)
    axes => axesi
  end subroutine set_axes

  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed

  double precision function norm2(v)
    double precision, intent(in), dimension(:) :: v
    norm2 = sum(v(:)*v(:))
  end function norm2

  double precision function norm(v)

    double precision, intent(in), dimension(:) :: v
    norm = dsqrt(norm2(v))
  end function norm

  integer function deltah(v,h2)
    double precision, intent(in), dimension(:) :: v
    double precision :: h2
    deltah = 0
    if (norm2(v) .le. h2) deltah = 1
  end function deltah

  double precision function varphi(v,h2)

    double precision, dimension(:) :: v
    double precision :: h2

    if (deltah(v,h2) .eq. 0) then
       varphi = 0
    else
       varphi = (1 - norm2(v)/h2)**3
    end if
  end function varphi

  double precision function phi(x,i,h2i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    double precision,intent(in) :: h2i
    double precision s,d,tmp
    double precision, dimension(size(x,1)) :: y
    integer j,ki

    y(:) = x(:)

    tmp = varphi(y,h2i)

    if (tmp .eq. 0) then
       phi = 0
       return
    end if

    d = 0
    do j=1,max_neighbors
       ki = node_neighbors2(i,j)
       if (ki .eq. 0) exit ! bad style!
       y(:) = x(:) + node_coordinates(i,:) - node_coordinates(ki,:) ! bad style!
       ! if (deltah(y,h2i) .eq. 0) cycle
       d = d + varphi(y,h2i)
    end do

    phi = 0
    if (d .ne. 0) phi = tmp / d
  end function phi

  double precision function testphijac()
    integer i,j
    double precision s,si, h2i
    double precision, dimension(3) :: a
    h2i = h ** 2

    s = 0.
    do j=1,numnodes
       call init_random_seed()
       call random_number(a)
       a(1)=2*(a(1)-0.5)
       a(2)=2*(a(2)-0.5)*dsqrt(1-a(1)**2)
       a(3)=sign(1.0D0,a(3)-0.5)*dsqrt(1-a(1)**2-a(2)**2)
       a(1)=a(1)*axes(1)
       a(2)=a(2)*axes(2)
       a(3)=a(3)*axes(3)
       si = 0.
       do i=1,numnodes
          si = si + phi(a(:)-node_coordinates(i,:),i,h2i)
       end do
       s = s + si
    end do
    testphijac = s
  end function testphijac

  subroutine normal_vector_stroke(numnodesi,neighbors)
    integer, intent(in),dimension(:,:) :: neighbors
    double precision, dimension(size(normal_coordinates,2)) :: tmpVector
    double precision, dimension(:), allocatable :: VSa
    integer, dimension(:), allocatable :: mask

    double precision VS, tmp, cand1, cand2

    integer i,numnodesi

    integer ind, count, icand, ncands

    tmpVector(:) = 0.

    do i=1, numnodesi
       ind = minloc(neighbors(i,:),dim=1) - 1
       allocate (mask(ind))
       mask = neighbors(i,1:ind)

       tmpVector(1) = sum(normal_coordinates(mask,1))
       tmpVector(2) = sum(normal_coordinates(mask,2))
       tmpVector(3) = sum(normal_coordinates(mask,3))

       nstroke_coordinates(i,:) = tmpVector(:) * (1.0 / norm(tmpVector))

       tmpVector(:) = nstroke_coordinates(i,:)

       ALLOCATE (VSa(ind))

       count = 0
       do while (count < 100)
          count = count + 1
          VSa(:) = tmpVector(1) * normal_coordinates(mask,1) + &
               tmpVector(2) * normal_coordinates(mask,2) + &
               tmpVector(3) * normal_coordinates(mask,3)
          VS = minval(VSa)

          icand = mask(minloc(VSa,dim = 1))
          ncands = ind
          cand2 = VS

          if (count .eq. 1) cand1 = cand2

          if (cand2 - cand1 .ge. -1.0E-17) then
             cand1 = cand2
             if (ncands .ne. 0) tmpVector(:) = tmpVector(:) + normal_coordinates(icand,:) / ncands
             tmp = 1.0 / norm(tmpVector) ! if 0, then bad logic here
             tmpVector(:) = tmpVector(:) * tmp
          end if
       end do
       nstroke_coordinates(i,:) = tmpVector(:)

       DEALLOCATE(VSa)
       DEALLOCATE(mask)
    end do
  end subroutine normal_vector_stroke

  double precision function dpair(i,j)
    integer, intent(in) :: i,j
    double precision, dimension(size(node_coordinates,2)) :: y

    y(:) = node_coordinates(i,:) - node_coordinates(j,:)

    dpair = norm(y)
    return
  end function dpair

  subroutine filter_neighbors(hvalue,neighbors,numnodesi)
    double precision, intent(in) :: hvalue
    integer, intent(in) :: numnodesi
    integer, intent(inout),dimension(:,:) :: neighbors

    integer i,j,ki
    double precision cand

    neighbors(:,:) = 0

    do i=1,numnodesi
       ki = 0
       do j=1,numnodesi
          cand = dpair(i,j)
          if ((i .eq. j) .or. (cand < hvalue)) then
             ki = ki + 1
             neighbors(i,ki) = j
          end if
       end do
    end do
  end subroutine filter_neighbors

end module modphi
