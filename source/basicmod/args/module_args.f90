module module_args
!-----------------------------------------------------------------------
! Module to manage command arguments
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 10/05/2013
!
! PRIVATE ATTRIBUTES:
!   args: array to save command arguments
!
! PUBLIC PROCEDURES:
!   get_arg: get i-th argument
!   is_arg: returns true when the argument is present
!   get_post_arg: returns the value after the required argument
!
! PRIVATE PROCEDURES:
!   get_args: get command arguments
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error
use module_convers, only: string
implicit none

!Private variables
character(maxpath), allocatable, private :: args(:)

!Private procedures
private :: get_args

contains

!-----------------------------------------------------------------------
! get_arg: get i-th argument
!-----------------------------------------------------------------------
function get_arg(i) result(res)
integer, intent(in) :: i
character(maxpath) :: res
integer :: nargs, length, st

nargs = command_argument_count()
if (i <= 0 .or. i > nargs) call error('(module_args/get_arg) required argument is out of range')
call get_command_argument(i, res, length, st)
if (st /= 0) call error('(module_args/get_arg) command line argument '//trim(string(i))//' not readable')
end function

!-----------------------------------------------------------------------
! is_arg: returns true when the argument is present
!-----------------------------------------------------------------------
function is_arg(argument) result(res)
character(*), intent(in) :: argument
logical ::  res
integer :: i

if (.not.allocated(args)) call get_args()
!do i = 1, size(args,1)
!print*,i,'.-', trim(args(i))
!end do
!print*,' '
res = .false.
do i = 1, size(args, 1)
  if (args(i) == argument) then
    res = .true.
    return
  end if
end do
end function

!-----------------------------------------------------------------------
! get_post_arg: returns the value after the required argument
!-----------------------------------------------------------------------
function get_post_arg(argument) result(res)
character(*), intent(in) :: argument
character(maxpath) :: res
integer :: i

if (.not. allocated(args)) call get_args()
do i = 1, size(args, 1) - 1
  if (args(i) == argument) then
    res = args(i+1)
    return
  end if
end do
call error('(module_args/get_post_arg) command line argument '//trim(string(argument))//' not found')
end function

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! get_args: get command arguments
!-----------------------------------------------------------------------
subroutine get_args()
character(maxpath) :: cad
integer :: nargs, length, st, res, i

nargs = command_argument_count()
allocate(args(nargs), stat = res, errmsg = cad)
!print*,'nargs', nargs
if (res /= 0) call error('(module_args/get_args) unable to allocate variable args: '//trim(cad))
do i = 1, nargs
  call get_command_argument(i, args(i), length, st)
  if (st /= 0) call error('(module_args/get_args) command line argument '//trim(string(i))//' not readable')
!print*,i,'.--', trim(args(i))
!print*,' '

end do
end subroutine


end module
