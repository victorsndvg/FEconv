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
subroutine load_mphtxt(mphtxtfile, nel, nnod, nver, DI_, LNN, LNV, LNE, LNF, nn, mm, nrc, nra, nrv, z, nsd)
character(len=*),         intent(in)  :: mphtxtfile
integer,                  intent(out) :: nel, nnod, nver, DI_, LNN, LNV, LNE, LNF
integer, allocatable,     intent(out) :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64),allocatable, intent(out) :: z(:,:)
type(mphtxt)                          :: u
type(pmh_mesh)                        :: mphtxt_m
type(mfm_mesh), allocatable           :: m(:)
integer                               :: d, DIM



  !inital settings
  call report_option('level', 'stdout')
  !process MPHTXT file
  print*, 'OPEN MPHTXT FILE:',trim(mphtxtfile)
  call open_mphtxt(u, mphtxtfile)
  print*, 'READING MPHTXT FILE:'
  print*, '-------------------------'
  call read_mphtxt(u, m, mphtxt_m, DIM)

end subroutine
end module
