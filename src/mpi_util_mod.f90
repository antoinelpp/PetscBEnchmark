module mpi_util_mod
  use MPI

  implicit none
  integer                                :: code, rang, Nprocs, comm2d, Nprocs_down
  integer,parameter                      :: ndims = 2
  integer,dimension(ndims)               :: dims, coords
  logical,dimension(ndims)               :: period
  integer,dimension(:,:),allocatable     :: topogrid
  integer,dimension(:),allocatable             :: lx,ly
  integer                                      :: sx, sy, ex, ey


contains

  !2D Topology creation
  subroutine topology_creation(xmax, ymax)
    implicit none
    integer :: xmax, ymax

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


end module mpi_util_mod
