module module_feed
!-----------------------------------------------------------------------
! Module to fill records in a file
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 27/01/2012
!
! PUBLIC PROCEDURES:
!   set_string_length: sets the maximum length of the internal string
!   set_ending: sets the ending for the internal string
!   feed: saves the argument in a internal string until it is filled; then
!     it is emptied in a file
!   empty: sends the content of the internal string to a file
!-----------------------------------------------------------------------
use module_convers
implicit none

!Variables
integer, private :: string_length  = 80  !the internal string will be full at this length
character(1024), private :: line   = ' ' !internal string
character(1024), private :: ending = ' ' !ending for the internal string
integer, private :: len_ending     = 0   !length of ending

contains

!-----------------------------------------------------------------------
! set_string_length: sets the maximum length of the internal string
!-----------------------------------------------------------------------
subroutine set_string_length(s)
integer, intent(in) :: s

string_length = min(s, 1024)
end subroutine

!-----------------------------------------------------------------------
! set_ending: sets the ending for the internal string
!-----------------------------------------------------------------------
subroutine set_ending(s)
character(*), intent(in) :: s

ending = s
len_ending = len_trim(ending)
end subroutine

!--------------------------------------------------------------------
! feed: saves the argument in a internal string until it is filled; then
! it is emptied in a file
!--------------------------------------------------------------------
subroutine feed(iu, arg)
integer,      intent(in) :: iu
character(*), intent(in) :: arg
character(len(arg)) :: adj_arg
      
adj_arg = adjustlt(arg)
if (string_length < len_trim(line) + 1 + len_trim(adj_arg) + len_ending) then
  write(iu, '(a)') trim(line)//trim(ending)
  line = trim(adj_arg)
elseif (len_trim(line) == 0) then
  line = trim(adj_arg)
else
  line = trim(line)//' '//trim(adj_arg)
end if
end subroutine

!--------------------------------------------------------------------
! empty: sends the content of the internal string to a file
!--------------------------------------------------------------------
subroutine empty(iu)
integer, intent(in) :: iu
  
write(iu, '(a)') trim(line)
line = ' '
end subroutine

end module
