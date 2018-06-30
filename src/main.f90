program main
  use mpi_util_mod
  use utilisis_mod
  use case_mod

  implicit none



  integer                                      :: xmax, ymax

  double PRECISION,dimension(:),allocatable    :: PHI, RHS

  CALL get_argument()

  xmax = xmax_arg
  ymax = ymax_arg

  CALL mpi_initial()
  CALL topology_creation(xmax, ymax)


  !print *, "KSP TYPE is:", trim(ksp_type), "  CP type is ;", TRIM(pc_type)

  ! print *,"range:",rang,"Nproc:", Nprocs,"boundaries", sx, sy, ex, ey

  ! print *,"range:",rang, lx, ly


  CALL petsc_initialize(MPI_COMM_WORLD, ymax+1, xmax+1, ly, lx, dims(2), dims(1), ksp_type, pc_type)


  call set_rho(sx,ex,sy,ey,xmax,ymax)

  BLOCK
    integer :: Number_of_my_cells
    Number_of_my_cells = ( ex - sx + 1 )*( ey - sy + 1 )
    allocate(RHS(Number_of_my_cells))
    allocate(PHI(Number_of_my_cells))
  END BLOCK

  call pack_tabToVec(tabrho, RHS)

  PHI(:) = 0.0

  CALL petsc_solver(RHS, PHI, log_view )

  call pack_vecToTab(tabphi, PHI)

  call calculate_error(sx,ex,sy,ey)

  CALL petsc_destroy()

  deallocate(PHI)
  deallocate(RHS)
  CALL mpi_final()


contains

  subroutine pack_tabToVec(tab, vec)
    implicit none
    integer :: indx, xx, yy, i, j

    double PRECISION,dimension(:), intent(out)            :: vec
    double PRECISION,dimension(sx:ex,sy:ey), intent(in)   :: tab

    do j = 1, ey - sy +1 !absolute y value
      do i = 1, ex - sx + 1 !absolute x value
        indx =  i + (ex - sx + 1)*(j-1)  !indx of the cell
        xx = (i + sx - 1) ! center of cell
        yy = (j + sy - 1) ! center of cell
        vec(indx) = tab(xx,yy)
      end do
    end do

  end subroutine pack_tabToVec

  subroutine pack_vecToTab(tab, vec)
    implicit none
    integer :: indx, xx, yy, i, j

    double PRECISION,dimension(:), intent(in)              :: vec
    double PRECISION,dimension(sx:ex,sy:ey), intent(out)   :: tab

    do j = 1, ey - sy +1 !absolute y value
      do i = 1, ex - sx + 1 !absolute x value
        indx =  i + (ex - sx + 1)*(j-1)  !indx of the cell
        xx = (i + sx - 1) ! center of cell
        yy = (j + sy - 1) ! center of cell
        tab(xx,yy) = vec(indx)
      end do
    end do

  end subroutine pack_vecToTab
  !*************************************************************************************************************


end program main
