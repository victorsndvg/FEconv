module module_dataset_2411_fcnv
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
use basicmod
use module_dataset_fcnv
!use module_mesh
use module_pmh_fcnv, only: piece
implicit none

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! read: read dataset 2411
!-----------------------------------------------------------------------
subroutine read_2411(iu, pc)

integer,        intent(in) :: iu   !unit number for unvfile
type(piece), intent(inout) :: pc   !piece
real(real64), allocatable, dimension(:) :: x
integer :: ios, Field1, F2, F3, F4
logical :: fit(2)

  call info('Reading mesh coordinates ...')

  pc%dim = 3
  if(.not. allocated(x)) allocate(x(pc%dim))

  do
    if (is_dataset_delimiter(iu, back=.true.)) exit
  ! Record 1
    read (unit=iu, fmt='(4I10)', iostat = ios) &
    Field1, & !node label
    F2, &     !export coordinate system number
    F3, &     !displacement coordinate system number
    F4        !color
    if (ios /= 0) call error('dataset_2411/read, #'//trim(string(ios)))
  ! Record 2
    x = 0._real64
    read (unit=iu, fmt='(3E25.16)', iostat = ios) &
    x !Fields 1-3    -- node coordinates in the part coordinate system
    if (ios /= 0) call error('dataset_2411/read, #'//trim(string(ios)))
  ! copy to mm
    fit = [.true., .false.]
    call set(2, pc%z, x(1:pc%dim), Field1, fit)
    pc%nnod = max(pc%nnod, Field1)
  end do

  call reduce(pc%z, pc%dim, pc%nnod)

end subroutine

end module
