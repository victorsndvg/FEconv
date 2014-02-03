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


function mphtxt_get_type(desc) result(res)

  character(len=*), intent(in)  :: desc
  integer                       :: res
  integer                       :: nnod, nver, lnn, lnv, lne, lnf

  nnod=0; nver=0; lnn=0; lnv=0; lne=0; lnf=0

  if(trim(desc) == 'vtx') then
    nnod=1; nver=1; lnn=1; lnv=1; lne=0; lnf=0
  elseif(trim(desc) == 'edg') then
    nnod=2; nver=2; lnn=2; lnv=2; lne=1; lnf=0
  elseif(trim(desc) == 'tri') then
    nnod=3; nver=3; lnn=3; lnv=3; lne=3; lnf=0
  elseif(trim(desc) == 'quad') then
    nnod=4; nver=4; lnn=4; lnv=4; lne=4; lnf=0
  elseif(trim(desc) == 'tet') then
    nnod=4; nver=4; lnn=4; lnv=4; lne=6; lnf=4
  elseif(trim(desc) == 'prism') then
    ! Prism FE not supported
  elseif(trim(desc) == 'hex') then
    nnod=8; nver=8; lnn=8; lnv=8; lne=12; lnf=6
  elseif(trim(desc) == 'edg2') then
    nnod=3; nver=2; lnn=3; lnv=2; lne=1; lnf=0
  elseif(trim(desc) == 'quad2') then
    nnod=8; nver=4; lnn=8; lnv=4; lne=4; lnf=0
  elseif(trim(desc) == 'tet2') then
    nnod=10; nver=4; lnn=10; lnv=4; lne=6; lnf=4
  elseif(trim(desc) == 'prism2') then
    ! Quadratic prism FE not supported
  elseif(trim(desc) == 'hex2') then
    nnod=16; nver=8; lnn=16; lnv=8; lne=12; lnf=6
  endif

  res = check_fe(nnod, nver, lnn, lnv, lne, lnf) 


end function

end module
