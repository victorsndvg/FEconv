module module_feed_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to send records to a file choosing the maximum length of the file lines and the continuation character.  
!!
!! Procedure `feed` adds information to a internal record, separating it from the previous one with a blank space;
!! when the internal record is filled, `feed` dumps it into the file.  
!!
!! Procedure `empty` dumps the remaining information in internal record into the file.  
!!  __Example:__  
!! `program test`  
!! `use basicmod, only: set_string_length, set_ending, feed, empty, string`  
!! `implicit none`  
!! `integer :: i, k`  
!! `character(10) :: s(100)`  
!! `s = [(string(2*i), i=1,100)]`  
!! `call set_string_length(132) !sets the maximum length of the file lines to 132`  
!! `call set_ending('\') !sets the continuation character to '\'`  
!! `open (10, file='data.dat')`  
!! `do k = 1, 100`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`call feed(10, s(k)) !adds information to the internal record and dumps it when it is filled`  
!! `end do`  
!! `call empty(10) !dumps the remaining of the internal record into the file`  
!! `close(10)`  
!! `end program`  
! 
! PUBLIC PROCEDURES:
!   set_string_length: sets the maximum length of a line
!   set_ending: sets the continuation character
!   feed: adds information to the current record
!   empty: save the current record in a file.
!-----------------------------------------------------------------------
use module_convers_bmod
implicit none

!Variables
integer, private :: string_length  = 80  !! Maximum length of a file line.
character(1024), private :: line   = ' ' !! Internal record.
character(1024), private :: ending = ' ' !! Continuation character.
integer, private :: len_ending     = 0   !! Length of the continuation character.

contains

!-----------------------------------------------------------------------
! set_string_length
!-----------------------------------------------------------------------
subroutine set_string_length(s)
!! Sets the maximum length of a file line; 80 characters by default. The maximum allowed value is 1024.
integer, intent(in) :: s !! Maximum length of the file line.

string_length = min(s, 1024)
end subroutine

!-----------------------------------------------------------------------
! set_ending
!-----------------------------------------------------------------------
subroutine set_ending(s)
!! Sets the continuation character; blank space by default. It can have a length bigger than one. 
character(*), intent(in) :: s !! Continuation character.

ending = s
len_ending = len_trim(ending)
end subroutine

!--------------------------------------------------------------------
! feed
!--------------------------------------------------------------------
subroutine feed(iu, arg)
!! Adds information to a internal record, separating it from the previous one with a blank space; 
!! when the internal record is filled, `feed` dumps it into the file.
integer,      intent(in) :: iu  !! File unit.
character(*), intent(in) :: arg !! Information to add to the internal record.
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
! empty
!--------------------------------------------------------------------
subroutine empty(iu)
!! Dumps the remaining information in internal record into the file. 
integer, intent(in) :: iu !! File unit.

write(iu, '(a)') trim(line)
line = ' '
end subroutine

end module
