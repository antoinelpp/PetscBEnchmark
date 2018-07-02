module case_mod
  use mpi_util_mod
  implicit none

  double PRECISION,dimension(:,:), allocatable :: global_phi, global_rho
  double PRECISION,dimension(:,:), allocatable :: global_phiExact, global_rhoExact

  double PRECISION,dimension(:,:),allocatable  :: tabphi, tabrho


contains

  subroutine set_rho(sx,ex,sy,ey,xmax,ymax)

    implicit none
    integer :: sx,ex,sy,ey,xmax,ymax

    allocate(tabphi(sx:ex,sy:ey))
    allocate(tabrho(sx:ex,sy:ey))

    allocate(global_phiExact(0:xmax, 0:ymax))
    allocate(global_rhoExact(0:xmax, 0:ymax))
    allocate(global_phi(0:xmax, 0:ymax))
    allocate(global_rho(0:xmax, 0:ymax))

    global_rhoExact(:,:) = 2.0d0
    global_rhoExact(:,0) = 0.0d0
    global_rhoExact(:,ymax) = 0.0d0

    global_phiExact(:,:) = 0.0d0
    tabrho(sx:ex,sy:ey) = global_rhoExact(sx:ex,sy:ey)

    BLOCK
      integer  :: i,j

      do j = 0, xmax
        do i=1,ymax-1

          global_phiExact(j,i) = dble(i - ymax/2)**2 - dble(ymax/2)**2
        end do

      end do

    ! if (rang == 0) then
    !   print *,global_phiExact(xmax/2,:)
    ! end if

    end block

  end subroutine set_rho

  subroutine calculate_error(sx,ex,sy,ey)
    implicit none
    integer :: i,j, code
    double PRECISION  :: error, error_red, norm
    integer :: sx,ex,sy,ey

    ! if (rang == 0) then
    !   print *,global_phiExact(ex/2,sy:ey) - tabphi(ex/2, sy:ey)
    ! end if

    norm = SUM(abs(global_phiExact(ex/2,sy:ey)))
    error = SUM( (global_phiExact(sx:ex,sy:ey) - tabphi(sx:ex,sy:ey))**2)

    call MPI_REDUCE(error,error_red,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm2d,code)

    if (rang == 0) then
      print *, "the total error is :", error_red
      print *, "the relativ error is :", error_red/norm
    end if

  end subroutine calculate_error

end module case_mod
