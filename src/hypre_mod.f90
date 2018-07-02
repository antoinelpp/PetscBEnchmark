module hypre_mod

  use mpi_util_mod

  implicit none
  integer                                 :: Number_of_my_cells
  real(dbleprc),allocatable               :: rhs(:),vals(:),sol(:)
  integer                                 :: lb_j,ub_j, lb_i, ub_i, hypre_xmax, hypre_ymax
  integer :: Nstensil = 3
  integer, dimension(2) :: hypre_period
  integer, dimension(1:4) ::  hypre_topo

contains

  !*******************************************************************************************************

  ! subroutine Set_Poisson_Boundary
  !
  !   implicit none
  !
  !   integer,dimension(:),allocatable             :: lx,ly
  !
  !   !for prob dev
  !   integer :: sx_wallprob, ex_wallprob
  !   logical :: wall_prob
  !   real(dbleprc) :: V_wallprob
  !
  !   allocate(lx(0:dims(1)-1))
  !   allocate(ly(0:dims(2)-1))
  !
  !   !CATHODE
  !   if(ey == ymax) then
  !     tabgrid(:,ey )%Phi = 0.d0
  !     tabgrid(:,ey+1 )%Phi = 0.d0 !ghostcell
  !   end if
  !
  !   !ANODE
  !   if(sy == 0) then
  !       if (abs(Omega) .LT. 1d-6) then
  !         tabgrid(:,sy)%Phi = q*dT**2*V0/(mi*dX*dX)
  !       else
  !         tabgrid(:,sy)%Phi = q*dT**2*V0/(mi*dX*dX)*SIN(Omega*(compt_loops)*dT)
  !       endif
  !       tabgrid(:,sy-1)%Phi = 0.d0 !ghostcell
  !   endif
  !
  !   if(.NOT. periodicx) then
  !     if(sx == 0) then
  !       tabgrid(0 ,:)%Phi= 0.d0
  !       tabgrid(-1,:)%Phi = 0.d0
  !     endif
  !     if(ex == xmax) then
  !       tabgrid(xmax  ,:)%Phi = 0.d0
  !       tabgrid(xmax+1,:)%Phi = 0.d0
  !     endif
  !   endif
  !
  ! end subroutine Set_Poisson_Boundary

  !*******************************************************************************************************
  !initializing stencil values
  subroutine set_stencil_hypre(xmax_glob, ymax_glob)
    implicit none
    integer, intent(in) :: xmax_glob, ymax_glob
    real(dbleprc)                              :: Ve, Vw, Vn, Vs
    integer                                    :: i,j,k,indx, xx,yy
    real(dbleprc),allocatable                  :: sig_up(:), sig_do(:)

    hypre_xmax = xmax_glob
    hypre_ymax = ymax_glob

    hypre_period = (/hypre_xmax+1, 0 /)

    !get the number of cells in each domain
    !Number_of_my_cells  = ( ex - sx + 1 )*( ey - sy + 1 )
    lb_j = 1
    ub_j = ey - sy + 1
    lb_i = 1
    ub_i = ex - sx + 1
    if(sy == 0)    lb_j = lb_j + 1     !bottom
    if(ey == hypre_ymax) ub_j = ub_j - 1     !top

    Number_of_my_cells = (ub_i - lb_i + 1)*(ub_j - lb_j + 1)

    !allocate the arrays
    if(.NOT. allocated(vals)) allocate( vals(Number_of_my_cells*Nstensil) )

    !pack the shit up
    do j = lb_j, ub_j
      do i = lb_i, ub_i
        indx =  i + 1 - lb_i +  (ub_i - lb_i + 1)*(j-lb_j)
          k = 1 + ( indx - 1 ) * Nstensil
          vals(k  ) = -4.0d0 ! C
          vals(k+1) =  1.0d0 ! E
          vals(k+2) =  1.0d0 ! N

          ! if(compt_loops == 0 .OR. (compt_loops == compt_restart .AND. compt_restart /= 1)) then
          !   k = 1 + ( indx - 1 ) * Nstensil
          !   vals(k+2) = 0.0d0
          ! end if
      enddo
    enddo

  end subroutine set_stencil_hypre

  !*******************************************************************************************************
  subroutine set_rhs_hypre(tabrho, tabphi)
    implicit none
    integer :: i,j,indx, xx, yy
    double PRECISION,dimension(sx:ex,sy:ey)  :: tabphi, tabrho



    if(.NOT. allocated(rhs) ) allocate( rhs (Number_of_my_cells) )
    if(.NOT. allocated(sol )) allocate( sol (Number_of_my_cells) )

    !pack the shit up
    do j = lb_j, ub_j
      do i = lb_i, ub_i
        indx =  i + 1 - lb_i +  (ub_i - lb_i + 1)*(j-lb_j)

        xx = (i + sx - 1) ! center of cell
        yy = (j + sy - 1) ! center of cell

        rhs(indx) = tabrho(xx,yy)
        sol(indx) = tabphi(xx,yy)
      enddo
    enddo


        !deal with the system boundaries
        !bottom
        if( sy == 0 ) then
          j = lb_j  ! on bottom
          do i = lb_i, ub_i
            indx =  i - lb_i + 1 + (ub_i - lb_i + 1)*(j-lb_j)
            xx = (i + sx - lb_i) ! center of cell
            yy = 0               ! interface
            ! South:
            rhs(indx) = rhs(indx)  - tabphi(xx,yy)
          enddo
        endif
        !top
        if( ey == hypre_ymax )then
          j = ub_j
          do i = lb_i, ub_i
            indx =  i - lb_i + 1 + (ub_i - lb_i + 1)*(j-lb_j)


            xx = (i + sx - lb_i) ! center of cell
            yy = hypre_ymax         ! interface
            ! North:
            rhs(indx) = rhs(indx) - tabphi(xx,yy)
          enddo
        endif

  end subroutine set_rhs_hypre

  !*******************************************************************************************************
  !Unpacking the results
  subroutine hypre_unpack(tabphi)
    implicit none

    integer                                  :: i,j,indx, xx,yy
    double PRECISION,dimension(sx:ex,sy:ey)  :: tabphi

    !lower and upper j bounds
    lb_j = 1
    ub_j = ey - sy + 1
    lb_i = 1
    ub_i = ex - sx + 1

    if(sy == 0)    lb_j = lb_j + 1     !bottom
    if(ey == hypre_ymax) ub_j = ub_j - 1     !top

    !print *,"shape",SHAPE(sol), ub_j- lb_j,ex - sx,(ub_j- lb_j+1)*(ex - sx+1)
    !if(rang==0) print *,"shape",SHAPE(sol), (ub_j- lb_j+1),&
    !  (ex - sx+1),(ub_j- lb_j+1)*(ex - sx+1)
    do j = lb_j, ub_j
      do i = lb_i, ub_i
        indx =  i - lb_i + 1 + (ub_i - lb_i + 1)*(j-lb_j)
        xx = (i + sx - 1) ! center of cell
        yy = (j + sy - 1) ! center of cell
        tabphi(xx,yy) = sol(indx)
        if(ISNAN( sol(indx))) then
          print *,"probleme hypre",sx,xx,ex,sy,yy,ey
        endif
      enddo
    enddo
    !if(rang==0) print *,"hypre unpacked"

  end subroutine hypre_unpack

  !*******************************************************************************************************
end module hypre_mod
