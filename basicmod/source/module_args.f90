module module_args_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to manage command arguments.
!
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
use module_compiler_dependant_bmod, only: real64
use module_os_dependant_bmod, only: maxpath
use module_report_bmod, only: error, info
use module_convers_bmod, only: string, word_count, word
implicit none

!Private variables
character(maxpath), allocatable, private :: args(:) !! Character array to save arguments.

!Private procedures
private :: set_args_sca, set_args_vec, set_args_from_command

!Interfaces
interface set_args
  !! Manually save arguments in a private variable.
  module procedure set_args_sca
  module procedure set_args_vec
end interface

contains

!-----------------------------------------------------------------------
! set_args_sca (set_args)
!-----------------------------------------------------------------------
subroutine set_args_sca(str)
!! Manually save arguments given as space-separated words in a character scalar.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set_args`  
!! `implicit none`  
!! `call set_args('-x -v -f filename)`  
!! `end program`  
character(*), intent(in) :: str !! Character scalar containing the arguments given as space-separated words.
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
! set_args_vec (set_args)
!-----------------------------------------------------------------------
subroutine set_args_vec(str)
!! Manually save arguments given as components in an character array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set_args`  
!! `implicit none`  
!! `character(4) :: str(3) = ['-x  ', '-f  ', 'fun ']`  
!! `call set_args(str)`  
!! `end program`  
character(*), allocatable :: str(:) !Character array containing arguments given as components.
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
! args_count
!-----------------------------------------------------------------------
function args_count() result(res)
!! Gets the number of arguments.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: args_count`  
!! `implicit none`  
!! `print*, 'Number of arguments given in the command line: ', args_count()`  
!! `end program`  
!!
!! @note The procedure reads arguments from command line if they were not previously saved.
integer :: res !! Number of saved arguments.

if (.not. allocated(args)) call set_args_from_command()
res = size(args, 1)
end function

!-----------------------------------------------------------------------
! get_arg
!-----------------------------------------------------------------------
function get_arg(i) result(res)
!! Retrieves `i`-th argument.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: get_arg`  
!! `implicit none`  
!! `print*, 'Argument n.1: ', trim(get_arg(1))`  
!! `end program`  
!!
!! @note The procedure reads arguments from command line if they were not previously saved.
integer, intent(in) :: i   !! Position of the argument to retrieve.
character(maxpath)  :: res !! Saved argument in the `i`-th position.
integer :: nargs

if (.not. allocated(args)) call set_args_from_command()
nargs = size(args, 1)
if (i <= 0 .or. i > nargs) call error('(module_args/get_arg) required argument '//trim(string(i))//' is out of range')
res = args(i)
end function

!-----------------------------------------------------------------------
! is_arg
!-----------------------------------------------------------------------
function is_arg(argument) result(res)
!! Checks whether a given argument is among the saved arguments.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_arg`  
!! `implicit none`  
!! `print*, 'Argument ''-d'' is given? ', is_arg('-d')`  
!! `end program`  
!!
!! @note The procedure reads arguments from command line if they were not previously saved.
character(*), intent(in) :: argument !! Argument to search.
logical                  ::  res     !! `.true.` when the argument is present; `.false.` otherwise.
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
! get_post_arg
!-----------------------------------------------------------------------
function get_post_arg(argument) result(res)
!! Retrieves the saved argument after the given one.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_arg, get_post_arg`  
!! `implicit none`  
!! `if (is_arg('-d')) then`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`print*, 'Argument following ''-d'': ', trim(get_post_arg('-d'))`  
!! `end if`  
!! `end program`  
!!
!! @note The procedure reads arguments from command line if they were not previously saved.
character(*), intent(in) :: argument !! Argument to search.
character(maxpath)       :: res      !! Saved argument after the given one.
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
! set_args_from_command
!-----------------------------------------------------------------------
subroutine set_args_from_command()
!! Sets arguments in args from command line.  
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
