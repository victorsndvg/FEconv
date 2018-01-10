module module_os_dependant
!-----------------------------------------------------------------------
! Module to define entities that depend on operating system
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 12/01/2012
!
! PUBLIC CONSTANTS:
!   maxpath: In the Windows API, the maximum length for a path;
!            also defined in Linux for compatibility issues
!
! PUBLIC PROCEDURES:
!   set_os: fill variable 'os' with the actual operating system
!     Warning: Recognizes Windows with gfortran as 'linux'
!   get_os: return variable 'os'
!   slash: returns the symbol for folder separation 
!-----------------------------------------------------------------------
use module_compiler_dependant, only: dir_exists
implicit none

!Constants
integer, parameter, public :: maxpath = 260

!Variables
character(maxpath), private :: os  = 'linux' !current operating system ('linux'|'windows')

contains

!-----------------------------------------------------------------------
! set_os: sets the actual value of OS
!-----------------------------------------------------------------------
subroutine set_os(arg)
character(*), optional, intent(in) :: arg
  
if (present(arg)) then
  os = arg
elseif (dir_exists('c:\')) then
  os = 'windows' 
else
  os = 'linux'
endif
end subroutine

!-----------------------------------------------------------------------
! get_os: gets the actual value of OS
!-----------------------------------------------------------------------
function get_os() result(res)
character(maxpath) :: res

res = os
end function

!-----------------------------------------------------------------------
! slash: returns the symbol for folder separation
!--------------------------------------------------------------------
function slash() result(res)
character(1) :: res

select case(trim(os))
case('windows')
  res = '\'
case default
  res = '/'
end select
end function

end module
