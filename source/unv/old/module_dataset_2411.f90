module module_dataset_2411
!-----------------------------------------------------------------------
! Module for dataset_2411 class
! Last update: 15/07/2008
!
! Name:   Nodes - Double Precision
! Record 1:        FORMAT(4I10)
!                  Field 1       -- node label
!                  Field 2       -- export coordinate system number
!                  Field 3       -- displacement coordinate system number
!                  Field 4       -- color
! Record 2:        FORMAT(1P3D25.16)
!                  Fields 1-3    -- node coordinates in the part coordinate
!                                   system
!
! Records 1 and 2 are repeated for each node in the model.
!
! Example:
!
!     -1
!   2411
!        121         1         1        11
!    5.0000000000000000D+00   1.0000000000000000D+00   0.0000000000000000D+00
!        122         1         1        11
!    6.0000000000000000D+00   1.0000000000000000D+00   0.0000000000000000D+00
!     -1
!-----------------------------------------------------------------------
use module_COMPILER_DEPENDANT, only: real64
use module_ALLOC
use module_dataset
use module_mesh_unv
implicit none

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! read: read dataset 2411
!-----------------------------------------------------------------------
subroutine read_2411(iu, m)

integer,    intent(in)    :: iu   !unit number for unvfile
type(mfm_mesh), intent(inout) :: m    !mesh
real(real64), dimension(3) :: x
integer :: ios, Field1, F2, F3, F4
logical :: fit(2)

do
  if (is_dataset_delimiter(iu, back=.true.)) return
! Record 1
  read (unit=iu, fmt='(4I10)', iostat = ios) &
  Field1, & !node label
  F2, &     !export coordinate system number
  F3, &     !displacement coordinate system number
  F4        !color
  if (ios /= 0) call error('dataset_2411/read, #'//trim(string(ios)))
! Record 2
  read (unit=iu, fmt='(3E25.16)', iostat = ios) &
  x !Fields 1-3    -- node coordinates in the part coordinate system
  if (ios /= 0) call error('dataset_2411/read, #'//trim(string(ios)))
! copy to mm
  fit = [.true., .false.]
  call set(2, m%xd, x(1:3), Field1, fit)
  m%nd = max(m%nd, Field1)
end do

end subroutine

end module
