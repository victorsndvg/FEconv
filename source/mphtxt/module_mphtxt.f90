module module_mphtxt

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

use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report
use module_convers
use module_manage_mphtxt
use module_mesh
use module_pmh
use module_fe_database_pmh

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
  call report_option('level', 'stdout')
 
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
