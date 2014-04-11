!            _   _ _
!      _   _| |_(_) |___
!     | | | | __| | / __|
!     | |_| | |_| | \__ \
!      \__,_|\__|_|_|___/
!

double precision function sigmaij(i,j)
  integer, intent(in) :: i,j
  sigmaij = dsqrt(sigma(i)**2 + sigma(j)**2)
end function sigmaij

!                      _             _
!      __   _____  ___| |_ ___  _ __| |__
!      \ \ / / _ \/ __| __/ _ \| '__| '_ \
!       \ V /  __/ (__| || (_) | |  | |_) |
!        \_/ \___|\___|\__\___/|_|  |_.__/
!

double complex function vectorb(i)
  use dbsym
  integer, intent(in) :: i

  vectorb = cdexp((0,1)*k_wave*node_coordinates(i,3))*intphi_over(i)

end function vectorb

double complex function vectorb2(i)
  use dbsym
  integer, intent(in) :: i

  vectorb2 = (0,1)*k_wave*normal_coordinates(i,3)*cdexp((0,1)*k_wave*node_coordinates(i,3))*intphi_over(i)

end function vectorb2

!       __  __       _        _          _
!      |  \/  | __ _| |_ _ __(_)_  __   / \
!      | |\/| |/ _` | __| '__| \ \/ /  / _ \
!      | |  | | (_| | |_| |  | |>  <  / ___ \
!      |_|  |_|\__,_|\__|_|  |_/_/\_\/_/   \_\
!

double complex function matrixA6(i,j)
  use omp_lib
  use dbsym
  integer, intent(in) :: i,j
  double precision :: sigm
  double precision, dimension(dim_3d) :: x,y
  integer :: jj, kj
  double complex :: tmp

  sigm = sigmaij(i,j)
  x = node_coordinates(i,:)
  y = node_coordinates(j,:)

  if (i .eq. j) then
     if (use_int_neighbors_p) then
        do jj=1, max_neighbors
           kj = node_neighbors2(i,jj)
           if (kj .eq. i) matrixA6 = intphi_over(i) * int_neighbors2(i,jj)/(4*PI)
        end do
     else
        matrixA6 = intphi_over(i) * (gauss(i,4) + gauss(i,7))
     end if
  else
     do jj=1, max_neighbors
        kj = node_neighbors2(i,jj)
        if (kj .eq. j) exit
        kj = 0
     end do
     if (kj > 0) then
        matrixA6 = intphi_over(i)*int_neighbors2(i,jj)/(4*PI)
     else
        matrixA6 = intphi_over(i)*foldingG(dim_intG,node_coordinates(i,:),j,k_wave)
      end if
  end if
end function matrixA6

double complex function matrixA(i,j)
  use omp_lib
  use dbsym
  integer, intent(in) :: i,j
  double precision :: sigm
  double precision, dimension(dim_3d) :: x,y

  sigm = sigmaij(i,j)
  x = node_coordinates(i,:)
  y = node_coordinates(j,:)

  counter = counter + 1
  if (i .eq. j) then
     matrixA = (intphi_over(i)**2)*intphi_under_micro(i)/(4*PI)/intphi_micro(i) !*limA(sigm,k_wave)
  else
     matrixA = intphi_over(i)*intphi_over(j)*Amn(x,y,k_wave)
  end if
end function matrixA

double complex function matrixA2(i,j)
  use dbsym
  integer, intent(in) :: i,j
  integer :: l
  double precision :: sigm, nstar, s
  double precision, dimension(3) :: x,y

  nstar = sum(normal_coordinates(i,:)*(node_coordinates(i,:)-node_coordinates(j,:))*intphi_over(j))

  sigm = sigmaij(i,j)
  x = node_coordinates(i,:)
  y = node_coordinates(j,:)

  if (i .eq. j) then
     matrixA2 = (gauss(j,1))*(intphi_over(i))
  else
     ! print *, Bmn(x,y,k_wave)
     matrixA2 = intphi_over(i)*nstar*Bmn(x,y,k_wave)
  end if
end function matrixA2

double complex function matrixA3(i,j)
  use omp_lib
  use dbsym
  integer, intent(in) :: i,j
  double precision :: sigm
  double precision, dimension(dim_3d) :: x,y

  sigm = sigmaij(i,j)
  x = node_coordinates(i,:)
  y = node_coordinates(j,:)

  if (i .eq. j) then
     matrixA3 = intphi_over(i)*(gauss(i,4) - gauss(i,3))
  else
     matrixA3 = intphi_over(i)*intphi_over(j)*Amn(x,y,k_wave)
  end if
end function matrixA3

double complex function matrixA_sigm(i,j)
  use omp_lib
  use dbsym

  integer, intent(in) :: i,j
  double precision :: sigm
  double precision, dimension(dim_3d) :: x,y

  sigm = sigmaij(i,j)
  x = node_coordinates(i,:)
  y = node_coordinates(j,:)

  if (i .eq. j) then
     matrixA_sigm = (intphi_over(i)**2)*limA(sigm,k_wave)
  else
     matrixA_sigm = intphi_over(i)*intphi_over(j)*A(x,y,sigm,k_wave)
  end if
end function matrixA_sigm

!                                        _                 _
!         __ _ _ __  _ __  _ __ _____  _(_)_ __ ___   __ _| |_ ___ _   _
!        / _` | '_ \| '_ \| '__/ _ \ \/ / | '_ ` _ \ / _` | __/ _ \ | | |
!       | (_| | |_) | |_) | | | (_) >  <| | | | | | | (_| | ||  __/ |_| |
!        \__,_| .__/| .__/|_|  \___/_/\_\_|_| |_| |_|\__,_|\__\___|\__,_|
!             |_|   |_|
!

double complex function approximateu4(x)
  use omp_lib
  use dbsym
  double precision, intent(in), dimension(:) :: x
  integer :: i
  double complex, dimension(:), allocatable :: s

  allocate(s(numnodes))

  call OMP_SET_NUM_THREADS(4)

  y_tmp(:) = x(:)
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(SHARED) PRIVATE(nt)
  do i=1,numnodes
     nt = OMP_GET_THREAD_NUM()
     if (norm(node_coordinates(i,:)-x) < 1.5*hval) then
     call integrate(                &
          fjacobian,                &
          nstroke_coordinates(i,:), &
          node_coordinates(i,:),    &
          i,                        &
          quadphi_over(:,1),        &
          quadphi_over(:,2),        &
          nt)
     s(i) = q_density(i) * folding(i,i,fAmn,dim_quad,nt)
     else
        s(i) = q_density(i) * foldingG(dim_intG,x,i,k_wave)
     end if
  end do
  !$OMP END PARALLEL DO

  approximateu4 = sum(s(:))

  deallocate(s)
end function approximateu4


double complex function approximateu(x)
  use omp_lib
  use dbsym
  double precision, intent(in), dimension(:) :: x
  integer :: i
  double complex, dimension(:), allocatable :: s

  allocate(s(numnodes))

  call OMP_SET_NUM_THREADS(4)

  !$OMP PARALLEL DO &
  !$OMP DEFAULT(SHARED) PRIVATE(v)
  do i=1,numnodes
     s(i) = intphi_over(i)*q_density(i)*Amn(x,node_coordinates(i,:),k_wave)
  end do
  !$OMP END PARALLEL DO

  approximateu = sum(s(:))

  deallocate(s)
end function approximateu

double complex function approximateu_sigm(x)
  use omp_lib
  use dbsym
  double precision, intent(in), dimension(:) :: x
  double precision :: sigm
  integer :: i
  double complex, dimension(:), allocatable :: s

  allocate(s(numnodes))

  call OMP_SET_NUM_THREADS(4)

  !$OMP PARALLEL DO &
  !$OMP DEFAULT(SHARED) PRIVATE(v)
  do i=1,numnodes
     sigm = sigma(i)
     s(i) = intphi_over(i)*q_density(i)*A(x,node_coordinates(i,:),sigm,k_wave)
  end do
  !$OMP END PARALLEL DO

  approximateu_sigm = sum(s(:))

  deallocate(s)
end function approximateu_sigm
