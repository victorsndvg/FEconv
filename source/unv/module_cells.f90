module module_cells_fcnv
!-----------------------------------------------------------------------
! Module for cell management
! Last update: 03/04/2010   
!
! cells is a matrix where:
!  - cells(1,i) is the dimension of the i-th cell in UNV file
!  - cell(2,i) is the index k in which the i-th cell is stored 
!    in the Modulef mesh of dimension cells(1,i)
!-----------------------------------------------------------------------
implicit none

!Variables
integer, public, allocatable, dimension(:,:) :: cells

end module
