module module_feconv
!-----------------------------------------------------------------------
! Module to convert between several mesh and FE field formats
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 10/05/2013
!
! PUBLIC PROCEDURES:
! convert: converts between several mesh and FE field formats
! is_arg: returns true when the argument is present
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error
use module_convers, only: adjustlt, lcase
use module_files, only: get_unit
use module_args, only: get_arg, is_arg, get_post_arg
use module_transform, only: lagr2l2, lagr2rt, lagr2nd
use module_cuthill_mckee, only: cuthill_mckee
use module_ansys, only: load_ansys
use module_unv, only: load_unv,save_unv
use module_patran, only: load_patran
use module_mfm, only: load_mfm, save_mfm
use module_mum, only: load_mum, save_mum
use module_vtu, only: save_vtu, type_cell
use module_mphtxt, only: load_mphtxt,save_mphtxt
use module_pf3, only: load_pf3,save_pf3
!use module_tra, only: load_tra,save_tra
use module_pmh
implicit none

!PMH structure
type(pmh_mesh) :: pmh

!Variables for MFM format
integer :: nel = 0 !global number of elements
integer :: nnod = 0 !global number of nodes
integer :: nver = 0 !global number of vertices
integer :: dim = 0 !space dimension
integer :: lnn = 0 !local number of nodes
integer :: lnv = 0 !local number of vertices
integer :: lne = 0 !local number of edges
integer :: lnf = 0 !local number of faces
integer, allocatable, dimension(:,:) :: nn !nodes index array
integer, allocatable, dimension(:,:) :: mm !vertices index array
integer, allocatable, dimension(:,:) :: nrv !vertices reference array
integer, allocatable, dimension(:,:) :: nra !edge reference array
integer, allocatable, dimension(:,:) :: nrc !face reference array
real(real64), allocatable, dimension(:,:) :: z !vertices coordinates array
integer, allocatable, dimension(:) :: nsd !subdomain index array

contains

!-----------------------------------------------------------------------
! convert: converts between several mesh and FE field formats
!-----------------------------------------------------------------------
subroutine convert()
character(maxpath) :: infile = ' ', inext = ' ', outfile = ' ', outext = ' '
integer :: p, nargs

!find infile and out file among arguments
nargs = command_argument_count()
infile = get_arg(nargs-1); p = index( infile, '.', back = .true.); inext = infile(p+1:len_trim( infile))
outfile = get_arg(nargs); p = index(outfile, '.', back = .true.); outext = outfile(p+1:len_trim(outfile))

!checks
if (len_trim(infile) == 0) call error('(module_feconv/fe_conv) unable to find input file.')
if (len_trim(outfile) == 0) call error('(module_feconv/fe_conv) unable to find output file.')
if (len_trim(inext) == 0) call error('(module_feconv/fe_conv) unable to find input file extension.')
if (len_trim(outext) == 0) call error('(module_feconv/fe_conv) unable to find output file extension.')
select case (trim(adjustlt(outext))) !check outfile extension now (avoid reading infile when outfile is invalid)
case('mfm', 'mum', 'vtu', 'mphtxt', 'unv')
  continue
case default
  call error('(module_feconv/fe_conv) output file extension not implemented: '//trim(adjustlt(outext)))
end select

!Change PMh mesh tolerance
if (is_arg('-t')) then
  pmh%ztol = dble(get_post_arg('-t'))
end if    

!read mesh
if (trim(adjustlt(inext)) /= 'unv' .and. is_arg('-is')) call error('(module_feconv/fe_conv) only UNV input files can '//&
&'manage -is option.')
select case (trim(lcase(adjustlt(inext))))
case('mfm')
  print '(a)', 'Loading MFM mesh file...'
  call load_mfm( infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('mum')
  print '(a)', 'Loading MUM mesh file...'
  call load_mum( infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('msh')
  print '(a)', 'Loading ANSYS mesh file...'
  call load_ansys(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('unv')
  print '(a)', 'Loading UNV mesh file...'
  call load_unv(infile, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, is_arg('-is'))
  print '(a)', 'Done!'
case('bdf')
  print '(a)', 'Loading MD Nastran input file...'
  call load_patran(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('mphtxt')
  print '(a)', 'Loading COMSOL mesh file...'
  call load_mphtxt(infile, pmh)
  print '(a)', 'Done!'
case('pf3')
  print '(a)', 'Loading FLUX mesh file...'
  call load_pf3(infile, pmh)
  print '(a)', 'Done!'
case default
  call error('(module_feconv/fe_conv) input file extension not implemented: '//trim(adjustlt(inext)))
end select
!print '(a)', 'FE type of the input mesh: '//trim(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))

!transform
if (is_arg('-l2')) then
print '(/a)', 'Converting Lagrange P1 mesh into Lagrange P2 mesh...'
  call lagr2l2(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
elseif (is_arg('-rt')) then
print '(/a)', 'Converting Lagrange P1 mesh into Raviart-Thomas (face) mesh...'
  call lagr2rt(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
elseif (is_arg('-nd')) then
print '(/a)', 'Converting Lagrange P1 mesh into Whitney (edge) mesh...'
  call lagr2nd(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
end if

!bandwidth optimization
if (is_arg('-cm')) then
call cuthill_mckee(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, z)
end if

!save mesh
select case (trim(adjustlt(outext)))
case('mfm')
  print '(/a)', 'Saving MFM mesh file...'
  if(allocated(pmh%pc)) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call save_mfm(outfile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('mum')
  print '(/a)', 'Saving MUM mesh file...'
  if(allocated(pmh%pc)) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call save_mum(outfile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('vtu')
  print '(/a)', 'Saving VTU mesh file...'
  if(allocated(pmh%pc)) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call save_vtu(outfile, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('mphtxt')
  print '(/a)', 'Saving COMSOL mesh file...'
  if(.not. allocated(pmh%pc)) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_mphtxt(outfile, pmh)
  print '(a)', 'Done!'
case('unv')
  print '(/a)', 'Saving I-DEAS UNV mesh file...'
  if(.not. allocated(pmh%pc)) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_unv(outfile, get_unit(), pmh)
  print '(a)', 'Done!'
end select !case default, already checked before reading infile
end subroutine

end module
