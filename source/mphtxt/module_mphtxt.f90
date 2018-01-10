module module_mphtxt_fcnv

!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! load_mphtxt: loads a mesh from a MPHTXT format file
! save_mphtxt: loads a mesh from a MPHTXT format file
!-----------------------------------------------------------------------

use basicmod
use module_manage_mphtxt_fcnv
!use module_mesh
use module_pmh_fcnv
use module_fe_database_pmh_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! load_mphtxt(filename, pmh): read a MPHTXT file
!-----------------------------------------------------------------------
! filename: file name
! pmh:      PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine load_mphtxt(filename, pmh)
  character(len=*), intent(in)          :: filename
  type(pmh_mesh), intent(inout)         :: pmh
  type(mphtxt)                          :: u
  integer                               :: dim

  ! Inital settings
  !call report_option('output', 'std')

  ! Open mphtxt file and reads the mesh
  call info('Reading mphtxt file ...')
  call open_mphtxt(u, filename)
  call read_mphtxt(u, pmh, dim)
  call close_mphtxt(u)

end subroutine

!-----------------------------------------------------------------------
! save_mphtxt(filename, pmh): write a MPHTXT file
!-----------------------------------------------------------------------
! filename: file name
! pmh:      PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine save_mphtxt(filename, pmh)
  character(len=*), intent(in) :: filename
  type(pmh_mesh), intent(inout) :: pmh
  type(mphtxt) :: u

  call open_mphtxt(u, filename,'unknown')
  call info('Writing mphtxt file ...')
  call write_mphtxt(u, pmh)
  call close_mphtxt(u)

end subroutine

end module
