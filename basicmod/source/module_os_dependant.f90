module module_os_dependant_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to define entities that depend on operating system.
!!
!! @note In the Windows API, the maximum length for a path is about 260 characters; we define variable `maxpath` to this value and 
!! we recommend to use it even in non-Windows operating systems for compatibility issues.  
!!
!! @warning Be aware that [[set_os(subroutine)]] compiled with compiler GFortran under MinGW classifies 
!! MS Windows as `'non-windows'`. This could be harmless since this compiler accepts linux-like paths.  
! 
! PUBLIC CONSTANTS:
!   maxpath: In the Windows API, the maximum length for a path; also defined in Linux for compatibility issues 
! 
! PUBLIC PROCEDURES:
!   set_os: fill variable 'os' with the actual operating system
!     Warning: Recognizes Windows with GFortran as 'non-windows'
!   get_os: return variable 'os'
!   slash: returns the symbol for folder separation
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: dir_exists, error_unit
implicit none

!Constants
integer, parameter, public :: maxpath = 260 !! Maximum path length (260 characters).

!Variables
character(maxpath), private :: os  = ' ' !! Actual operating system (`windows`or `non-windows`). By default, ` `.

contains

!-----------------------------------------------------------------------
! set_os
!-----------------------------------------------------------------------
subroutine set_os(arg)
!! Inquire the actual operating system.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set_os`  
!! `implicit none`  
!! `call set_os()`  
!! `end program`  
!!
!! @note When the argument is missing, this procedure searches folder `'c:\'` to set the operating system.  
!!
!! @warning Be aware that [[set_os(subroutine)]] compiled with compiler GFortran under MinGW classifies 
!! MS Windows as `'non-windows'`. This could be harmless since this compiler accepts linux-like paths.  
character(*), optional, intent(in) :: arg !! Operating system; when it is given, valid values are `windows` and `non-windows`.
integer :: ios

if (present(arg)) then
  select case(arg)
  case('windows', 'non-windows')
    os = arg
  case default
    write(error_unit,'(a)') 'ERROR: (module_os_dependant/set_os) Only valid choices for operating system are &
    &''windows'' and ''non-windows'''
    stop 1
  end select
elseif (dir_exists('c:\', ios)) then
  os = 'windows'
else
  os = 'non-windows'
endif
end subroutine

!-----------------------------------------------------------------------
! get_os: gets the actual value of OS
!-----------------------------------------------------------------------
function get_os() result(res)
!! Returns the actual operating system.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: get_os`  
!! `implicit none`  
!! `print*, get_os()`  
!! `end program`  
!!
!! @note If the operating system has not been saved, procedure [[set_os(subroutine)]] is called.  
character(maxpath) :: res !! Operating system: `windows` or `non-windows`.

if (os == ' ') call set_os()
res = os
end function

!-----------------------------------------------------------------------
! slash
!--------------------------------------------------------------------
function slash() result(res)
!! Returns the symbol for folder separation in the actual operating system.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set_os, slash`  
!! `implicit none`  
!! `call set_os()`  
!! `print*, slash()`  
!! `end program`  
!!
!! @note Call `set_os` to inquire the actual operating system.  
character(1) :: res !! Symbol for folder separation: `\` for MS Windows and `/` otherwise.

select case(trim(get_os()))
case('windows')
  res = '\'
case default
  res = '/'
end select
end function

end module
