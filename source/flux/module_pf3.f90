module module_pf3_fcnv

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! load_pf3: loads a mesh from a PF3 format file
! save_pf3: loads a mesh from a PF3 format file
!-----------------------------------------------------------------------

use basicmod
use module_manage_pf3_fcnv
!use module_mesh
use module_pmh_fcnv
use module_fe_database_pmh_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! load_pf3(filename, pmh): read a PF3 file
!-----------------------------------------------------------------------
! filename: file name
! pmh:      PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine load_pf3(filename, pmh)
  character(len=*), intent(in)          :: filename
  type(pmh_mesh), intent(inout)         :: pmh
  type(pf3)                             :: u
  integer                               :: dim

  ! Inital settings
  !call report_option('level', 'stdout')

  ! Open PF3 file and reads the mesh
  call info('Reading PF3 file ...')
  call open_pf3(u, filename)
  call read_pf3(u, pmh, dim)
  call close_pf3(u)

end subroutine

!-----------------------------------------------------------------------
! save_pf3(filename, pmh): write a PF3 file
!-----------------------------------------------------------------------
! filename: file name
! pmh:      PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine save_pf3(filename, pmh, infield, outfield, path, param)
  character(len=*),          intent(in) :: filename
  type(pmh_mesh),         intent(inout) :: pmh
  character(*), allocatable, intent(in) :: infield(:)  ! In field names
  character(*), allocatable, intent(in) :: outfield(:) ! Out field names
  character(*),              intent(in) :: path !file names
  real(real64), optional,    intent(in) :: param
  type(pf3) :: u

  call open_pf3(u, filename,'unknown')
  call info('Writing PF3 file ...')
  call write_pf3(u, pmh, infield, outfield, path, param)
  call close_pf3(u)

end subroutine

end module
