module module_mphtxt
!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 08/01/2014
!
! PUBLIC PROCEDURES:
!   load_mphtxt: loads a mesh from a MPHTXT format file
!   NOTES: 1) mesh must be composed of a single FE type
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
! load_mphtxt: read a MPHTXT file
!-----------------------------------------------------------------------
subroutine load_mphtxt(filename, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  character(len=*),         intent(in)   :: filename
  integer,                  intent(out)  :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
  integer, allocatable,     intent(out)  :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
  real(real64),allocatable, intent(out)  :: z(:,:)
  type(pmh_mesh),           intent(inout):: pmh
  type(mphtxt)                           :: u
  type(mfm_mesh), allocatable            :: m(:)


  ! Inital settings
  call report_option('level', 'stdout')
 
  ! Open mphtxt file and reads the mesh
  call info('Reading mphtxt file ...')
  call open_mphtxt(u, filename)
  call read_mphtxt(u, m, pmh, dim)
  call close_mphtxt(u)

  ! Build vertex connectivities and coordinates from node information
  call build_vertices(pmh)

end subroutine

!-----------------------------------------------------------------------
! save_mphtxt: write a MPHTXT file
!-----------------------------------------------------------------------
subroutine save_mphtxt(filename, pmh)
  character(len=*),         intent(in)     :: filename
  type(pmh_mesh),           intent(inout)  :: pmh
  type(mphtxt)                             :: u

  call open_mphtxt(u, filename,'unknown')
  call write_mphtxt(u, pmh)
  call close_mphtxt(u)

end subroutine

end module
