module utilisis_mod

  implicit none

  character(len=10) :: ksp_type, pc_type
  integer :: xmax_arg, ymax_arg, log_view=0




contains

  subroutine get_argument
  !this function reads all of the argument of the command

  implicit none

  !restart value
  character(len=100)                :: arg,value
  integer                           :: i,Narg, ios

  ! DEFAULT arguments
  pc_type = "jacobi"
  ksp_type = "cg"
  xmax_arg = 100
  ymax_arg = 200

  Narg = command_argument_count()
  !print *,Narg
  if(modulo(Narg,2) /= 0) stop "Argument must be even (-<opt> <valeur>)"

  do i=1,Narg,+2
    !print *,'i',i
    !getting the arguments in the execution command line
    call get_command_argument(i,arg)
    call get_command_argument(i+1,value)
    if (LEN_TRIM(arg) == 0) stop "erreur argument of command"

    select case (trim(arg))
    case('-ksp_type')
      ksp_type = trim(value)

    case('-pc_type')
      ! print '(a)', '  -i, --input_dir   set the value of the input directory'
      pc_type = trim(value)

    case('-log_view')
      ! print '(a)', '  -i, --input_dir   set the value of the input directory'
      read( value, '(I3)', iostat=ios) log_view

    case('-xmax')
      ! print '(a)', '  -i, --input_dir   set the value of the input directory'
      read( value, '(I10)', iostat=ios) xmax_arg

    case('-ymax')
      ! print '(a)', '  -i, --input_dir   set the value of the input directory'
      read( value, '(I10)', iostat=ios) ymax_arg

    case default
      stop "Argument de commande invalide"// trim(arg)
    end select
  end do


end subroutine get_argument

end module utilisis_mod
