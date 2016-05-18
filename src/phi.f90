module modphi
  use params
contains
  
  include 'set_params.f90'
  
  double precision function norm2(v)
    double precision, intent(in), dimension(:) :: v
    norm2 = sum(v(:)*v(:))
  end function norm2

  double precision function norm(v)

    double precision, intent(in), dimension(:) :: v
    norm = dsqrt(norm2(v))
  end function norm

  integer function deltah(v,hval2)
    double precision, intent(in), dimension(:) :: v
    double precision :: hval2
    deltah = 0
    if (norm2(v) .le. hval2) deltah = 1
  end function deltah

  double precision function varphi(v,hval2)

    double precision, dimension(:) :: v
    double precision :: hval2

    if (deltah(v,hval2) .eq. 0) then
       varphi = 0
    else
       varphi = (1 - norm2(v)/hval2)**3
    end if
  end function varphi

  double precision function phi(x,i,hval2i)
    integer, intent(in) :: i
    double precision, intent(in), dimension(:) :: x
    double precision,intent(in) :: hval2i
    double precision s,d,tmp
    double precision, dimension(size(x,1)) :: y
    integer j,ki

    y(:) = x(:)

    tmp = varphi(y,hval2i)

    if (tmp .eq. 0) then
       phi = 0
       return
    end if

    d = 0
    do j=1,max_neighbors
       ki = node_neighbors2(i,j)
       if (ki .eq. 0) exit ! bad style!
       y(:) = x(:) + node_coordinates(i,:) - node_coordinates(ki,:) ! bad style!
       ! if (deltah(y,hval2i) .eq. 0) cycle
       d = d + varphi(y,hval2i)
    end do

    phi = 0
    if (d .ne. 0) phi = tmp / d
  end function phi

  double precision function test_phi(a)
    double precision, intent(in), target :: a(:,:)
    integer i,j
    double precision s,si, hval2i
    hval2i = hval ** 2

    s = 0.
    do j=1,numnodes
       si = 0.
       do i=1,numnodes
          si = si + phi(a(j,:)-node_coordinates(i,:),i,hval2i)
       end do
       s = s + si
    end do
    test_phi = s
  end function test_phi

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
