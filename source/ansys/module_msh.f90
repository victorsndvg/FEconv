module module_msh_fcnv

!-----------------------------------------------------------------------
! Module to manage MSH (Ansys) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! load_msh: loads a mesh from a MSH format file
! save_msh: loads a mesh from a MSH format file
!-----------------------------------------------------------------------

use basicmod
!use module_mesh
use module_pmh_fcnv
use module_fe_database_pmh_fcnv
use module_manage_msh_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! load_msh(filename, pmh): read a MSH file
!-----------------------------------------------------------------------
! filename: file name
! pmh:      PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine load_msh(filename, pmh)
  character(len=*), intent(in)          :: filename
  type(pmh_mesh), intent(inout)         :: pmh
  type(msh)                             :: u

  ! Inital settings
  !call report_option('level', 'stdout')

  ! Open msh file and reads the mesh
  call info('Reading Ansys MSH file ...')
  call open_msh(u, filename)
  call read_msh(u, pmh)
  call close_msh(u)

end subroutine

!-----------------------------------------------------------------------
! save_msh(filename, pmh): write a MSH file
!-----------------------------------------------------------------------
! filename: file name
! pmh:      PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine save_msh(filename, pmh)
  character(len=*), intent(in) :: filename
  type(pmh_mesh), intent(inout) :: pmh
  type(msh) :: u

  call open_msh(u, filename,'unknown')
  call info('Writing Ansys MSH file ...')
  call write_msh(u, pmh)
  call close_msh(u)

end subroutine

end module
