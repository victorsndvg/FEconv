module module_dataset
!-----------------------------------------------------------------------
! Module for datasets 
!
! In a previous version, dataset was an abstract type. We changed this 
! since it appears that was incompatible with the allocatable member feL
!
! Last update: 04/04/2010   
!-----------------------------------------------------------------------
use module_REPORT
use module_CONVERS
implicit none

contains

!***********************************************************************
! MODULE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! is_dataset_delimiter: check whether or not there is a dataset delimiter
!-----------------------------------------------------------------------
function is_dataset_delimiter(iu, back) result(res)

integer, intent(in) :: iu !unit number for unvfile
logical, intent(in), optional :: back
logical :: res
integer :: val, ios

res = .false.
read (unit=iu, fmt='(I10)', iostat=ios) val
if (ios /= 0) call error('dataset/is_dataset_delimiter/read, #'//trim(string(ios)))
!delimiter found
if (val == -1) then
  res = .true.
  return
end if
!backspace
if (present(back)) then; if (back) then
  backspace(unit=iu, iostat=ios)
  if (ios /= 0) call error('dataset/is_dataset_delimiter/backspace, #'//trim(string(ios)))
end if; end if
end function

end module
