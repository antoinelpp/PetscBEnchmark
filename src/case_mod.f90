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

    global_rhoExact(:,:) = 1.0d0

    BLOCK
      integer  :: i
      do i=0,ymax
        global_phiExact(:,i) = dble(i - xmax/2)**2 - dble(xmax/2)**2
      end do


    end block

  end subroutine set_rho


  subroutine calculate_error(sx,ex,sy,ey)
    implicit none
    integer :: i,j, code
    double PRECISION  :: error, error_red
    integer :: sx,ex,sy,ey


    error = SUM( global_phiExact(sx:ex,sy:ey) - tabphi(sx:ex,sy:ey))

    call MPI_REDUCE(error_red,error,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm2d,code)

    if (rang == 0) then
      print *, "the total error is :", error_red
    end if


  end subroutine calculate_error

end module case_mod
