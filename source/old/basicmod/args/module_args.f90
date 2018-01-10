module module_args
!-----------------------------------------------------------------------
! Module to manage command arguments
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 15/05/2014
!
! PRIVATE ATTRIBUTES:
!   args: array to save command arguments
!
! PUBLIC PROCEDURES:
!   set_args: manually set arguments in args (see also set_args_from_command)
!   args_count: get the number of arguments
!   get_arg: get i-th argument
!   is_arg: returns true when the argument is present
!   get_post_arg: returns the value after the required argument
!
! PRIVATE PROCEDURES:
!   set_args_from_command: set arguments in args from command line
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_convers, only: string, word_count, word
implicit none

!Private variables
character(maxpath), allocatable, private :: args(:)

!Private procedures
private :: set_args_sca, set_args_vec, set_args_from_command

!Interfaces
interface set_args; module procedure set_args_sca; end interface
interface set_args; module procedure set_args_vec; end interface

contains

!-----------------------------------------------------------------------
! set_args: manually set arguments in args
!-----------------------------------------------------------------------
subroutine set_args_sca(str)
character(*) :: str
integer :: nargs, res, i
character(maxpath) :: cad

if (allocated(args)) then
  call info('(module_args/set_args_sca) arguments args already allocated; deallocating first...')
  deallocate(args, stat = res, errmsg = cad)
  if (res /= 0) call error('(module_args/set_args_sca) Unable to deallocate args')
end if
nargs = word_count(str)
allocate(args(nargs), stat = res, errmsg = cad)
if (res /= 0) call error('(module_args/set_args_sca) unable to allocate variable args')
do i = 1, nargs
  args(i) = word(str, i)
end do
end subroutine

!-----------------------------------------------------------------------
! set_args: manually set arguments in args
!-----------------------------------------------------------------------
subroutine set_args_vec(str)
character(*), allocatable :: str(:)
integer :: nargs, res, i
character(maxpath) :: cad

if (.not. allocated(str)) call error('(module_args/set_args_vec) input argument str is not allocated')
if (allocated(args)) then
  call info('(module_args/set_args_vec) arguments args already allocated; deallocating first...')
  deallocate(args, stat = res, errmsg = cad)
  if (res /= 0) call error('(module_args/set_args_vec) Unable to deallocate args')
end if
nargs = size(str, 1)
allocate(args(nargs), stat = res, errmsg = cad)
if (res /= 0) call error('(module_args/set_args_vec) unable to allocate variable args')
do i = 1, nargs
  args(i) = str(i)
end do
end subroutine

!-----------------------------------------------------------------------
! args_count: get the number of arguments
!-----------------------------------------------------------------------
function args_count() result(res)
integer :: res

if (.not. allocated(args)) call set_args_from_command()
res = size(args, 1)
end function

!-----------------------------------------------------------------------
! get_arg: get i-th argument
!-----------------------------------------------------------------------
function get_arg(i) result(res)
integer, intent(in) :: i
character(maxpath) :: res
integer :: nargs

if (.not. allocated(args)) call set_args_from_command()
nargs = size(args, 1)
if (i <= 0 .or. i > nargs) call error('(module_args/get_arg) required argument '//trim(string(i))//' is out of range')
res = args(i)
end function

!-----------------------------------------------------------------------
! is_arg: returns true when the argument is present
!-----------------------------------------------------------------------
function is_arg(argument) result(res)
character(*), intent(in) :: argument
logical ::  res
integer :: i

if (.not.allocated(args)) call set_args_from_command()
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

if (.not. allocated(args)) call set_args_from_command()
do i = 1, size(args, 1) - 1
  if (args(i) == argument) then
    res = args(i+1)
    return
  end if
end do
call error('(module_args/get_post_arg) argument '//trim(string(argument))//' not found')
end function

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! set_args_from_command: set arguments in args from command line
!-----------------------------------------------------------------------
subroutine set_args_from_command()
character(maxpath) :: cad
integer :: nargs, length, st, res, i

nargs = command_argument_count()
allocate(args(nargs), stat = res, errmsg = cad)
if (res /= 0) call error('(module_args/set_args_from_command) unable to allocate variable args: '//trim(cad))
do i = 1, nargs
  call get_command_argument(i, args(i), length, st)
  if (st /= 0) call error('(module_args/set_args_from_command) command line argument '//trim(string(i))//' not readable')
end do
end subroutine

end module
