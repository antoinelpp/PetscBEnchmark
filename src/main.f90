program main
  use MPI
  use utilisis_mod

  implicit none

  integer                                :: code, rang, Nprocs, comm2d, Nprocs_down
  integer,parameter                      :: ndims = 2
  integer,dimension(ndims)               :: dims, coords
  logical,dimension(ndims)               :: period
  integer,dimension(:,:),allocatable     :: topogrid

  integer                                      :: xmax, ymax, sx, sy, ex, ey
  integer,dimension(:),allocatable             :: lx,ly
  double PRECISION,dimension(:),allocatable    :: PHI, RHS
  double PRECISION,dimension(:,:),allocatable  :: tabphi, tabrho

  CALL get_argument()

  xmax = xmax_arg
  ymax = ymax_arg

  CALL mpi_initial()
  CALL topology_creation()


  !print *, "KSP TYPE is:", trim(ksp_type), "  CP type is ;", TRIM(pc_type)

  ! print *,"range:",rang,"Nproc:", Nprocs,"boundaries", sx, sy, ex, ey

  ! print *,"range:",rang, lx, ly


  CALL petsc_initialize(MPI_COMM_WORLD, ymax+1, xmax+1, ly, lx, dims(2), dims(1), ksp_type, pc_type)

  allocate(tabphi(sx:ex,sy:ey))
  allocate(tabrho(sx:ex,sy:ey))

  call RANDOM_NUMBER(tabrho)
  tabrho(:,:) = 0

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
  !2D Topology creation
  subroutine topology_creation
    implicit none

    !no reorganization made by MPI
    logical,parameter          :: reOrganization = .false.
    logical                    :: condition
    integer,dimension(:),allocatable  :: temp_topogrid
    integer, dimension(4) :: valeurs
    !knowing the number of processes along each axes in function of the total number of processes
    dims(:) = (/ 0,0 /)

    call MPI_DIMS_CREATE(Nprocs,ndims,dims,code)

    !creation of the topology
    period(:) = (/  .TRUE., .FALSE. /)
    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, period, reOrganization, comm2d, code)

    !Find local coords
    call MPI_CART_COORDS(comm2d, rang, ndims, coords, code)

    !print *,"dims = ",rang, dims, "coords =", coords


    if(coords(1) == 0) then
      sx = 0
    else
      sx = (coords(1)*xmax)/dims(1)+1
    end if
    ex = ((coords(1)+1)*xmax)/dims(1)

    if(coords(2) == 0) then
      sy = 0
    else
      sy = (coords(2)*ymax)/dims(2)+1
    end if
    ey = ((coords(2)+1)*ymax)/dims(2)

    allocate(temp_topogrid(Nprocs*4))
    if(.NOT. allocated(topogrid)) allocate(topogrid(0:Nprocs-1,4))

    temp_topogrid(:) = 0
    topogrid(:,:)    = 0

    valeurs(1) = sx
    valeurs(2) = sy
    valeurs(3) = ex
    valeurs(4) = ey

    !Gathering the global coordinates on every processor. So each one know where every one is.
    call MPI_ALLGATHER(valeurs,4,MPI_INTEGER,temp_topogrid,4,MPI_INTEGER,comm2d,code)
    if(code /= 0) print*,'bhytthe dans le ALLGATHER de la topo'


    allocate(lx(0:dims(1)-1))
    allocate(ly(0:dims(2)-1))

    BLOCK
      integer :: k, j, i
      k=0; j=0
      do i=0,Nprocs-1
        do j=1,4
          topogrid(i,j) = temp_topogrid(4*i+j)
        end do
      end do

      k=0; j=0
      do i=0,Nprocs-1

        if ( topogrid(i,2)==0 ) then ! sy =0
          lx(k)=topogrid(i,3)-topogrid(i,1) + 1
          k=k+1
        end if

        if(topogrid(i,1)==0 ) then ! sx = 0
          ly(j)=topogrid(i,4)-topogrid(i,2) + 1
          j=j+1
        end if

      end do

    END BLOCK
    deallocate(temp_topogrid)


  end subroutine topology_creation


  !*************************************************************************************************************
  !MPI process initialization
  subroutine mpi_initial
    implicit none

    !MPI's initialization
    call MPI_INIT(code)

    !Getting the rank of each process
    call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

    !Number of processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,code)

  end subroutine mpi_initial

  !*************************************************************************************************************
  !Finalize MPI process
  subroutine mpi_final
    implicit none

    call MPI_FINALIZE(code)

  end subroutine mpi_final

end program main
