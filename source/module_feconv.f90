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
use module_alloc, only: set
use module_args, only: get_arg, is_arg, get_post_arg, args_count
use module_transform, only: lagr2l2, lagr2rt, lagr2nd, to_l1
use module_cuthill_mckee, only: cuthill_mckee
use module_msh, only: load_msh,save_msh
use module_unv, only: load_unv,save_unv
use module_patran, only: load_patran
use module_mfm, only: load_mfm, save_mfm
use module_mum, only: load_mum, save_mum
use module_vtu, only: save_vtu2, load_vtu, type_cell
use module_mphtxt, only: load_mphtxt,save_mphtxt
use module_pf3, only: load_pf3,save_pf3
!use module_tra, only: load_tra,save_tra
use module_field_database, only: FLDB, id_mesh_ext
use module_freefem, only: save_freefem
use module_pmh
implicit none

!PMH structure
type(pmh_mesh) :: pmh

!Variables for MFM format
integer :: nel  = 0 !global number of elements
integer :: nnod = 0 !global number of nodes
integer :: nver = 0 !global number of vertices
integer :: dim  = 0 !space dimension
integer :: lnn  = 0 !local number of nodes
integer :: lnv  = 0 !local number of vertices
integer :: lne  = 0 !local number of edges
integer :: lnf  = 0 !local number of faces
integer, allocatable, dimension(:,:) :: nn  !nodes index array
integer, allocatable, dimension(:,:) :: mm  !vertices index array
integer, allocatable, dimension(:,:) :: nrv !vertices reference array
integer, allocatable, dimension(:,:) :: nra !edge reference array
integer, allocatable, dimension(:,:) :: nrc !face reference array
real(real64), allocatable, dimension(:,:) :: z !vertices coordinates array
integer, allocatable, dimension(:) :: nsd !subdomain index array

logical :: is_pmh !true if the working mesh is PMH, false if is MFM

contains

!-----------------------------------------------------------------------
! convert: converts between several mesh and FE field formats
!-----------------------------------------------------------------------
subroutine convert()
character(maxpath) :: infile=' ', inmesh=' ', inext=' ', outfile=' ', outmesh=' ', outext=' '
character(maxpath), allocatable :: infield(:), outfield(:)
integer :: p, nargs
logical :: there_is_field

!find infile and outfile at the end of the arguments
nargs = args_count()
 infile = get_arg(nargs-1); p = index( infile, '.', back=.true.);  inmesh =  infile(1:p-1);  inext =  infile(p+1:len_trim( infile))
outfile = get_arg(nargs);   p = index(outfile, '.', back=.true.); outmesh = outfile(1:p-1); outext = outfile(p+1:len_trim(outfile))

!check mesh names and extensions
if (len_trim(infile)  == 0) call error('(module_feconv/fe_conv) unable to find input file.')
if (len_trim(outfile) == 0) call error('(module_feconv/fe_conv) unable to find output file.')
if (len_trim(inext)   == 0) call error('(module_feconv/fe_conv) unable to find input file extension.')
if (len_trim(outext)  == 0) call error('(module_feconv/fe_conv) unable to find output file extension.')
select case (trim(adjustlt(outext))) !check outfile extension now (avoid reading infile when outfile is invalid)
case('mfm', 'mum', 'vtu', 'mphtxt', 'unv', 'pf3', 'msh', 'mesh')
  continue
case default
  call error('(module_feconv/fe_conv) output file extension not implemented: '//trim(adjustlt(outext)))
end select
!check isoparametric option, for UNV only
if (trim(adjustlt(inext)) /= 'unv' .and. is_arg('-is')) call error('(module_feconv/fe_conv) only UNV input files can '//&
&'manage -is option.')
!options for mesh transformation (-l1, -l2, -rt, -nd and -cm) are incompatible with fields (-if, -of)
if ( (is_arg('-l1') .or. is_arg('-l2') .or. is_arg('-rt') .or. is_arg('-nd') .or. is_arg('-cm')) .and. &
     (is_arg('-if') .or. is_arg('-of')) ) call error('(module_feconv/fe_conv) options for mesh transformation (-l1, -l2, '//&
     &'-rt, -nd and -cm) are incompatible with fields (-if, -of).')
!set PMH mesh tolerance (all load procedures must consider intent(inout) for PMH argument)
if (is_arg('-t')) then 
  pmh%ztol = dble(get_post_arg('-t'))
end if    

!field selection
there_is_field = .true.
if (FLDB(id_mesh_ext(inext))%is_field_outside) then
  if (FLDB(id_mesh_ext(outext))%is_field_outside) then
    !infield and outfield are both mesh external
    if (is_arg('-if') .and. is_arg('-of')) then
      !there is -if, there is -of
      call set( infield, get_post_arg('-if'), 1, fit=.true.)
      call set(outfield, get_post_arg('-of'), 1, fit=.true.)
    elseif (.not. is_arg('-if') .and. .not. is_arg('-of')) then
      !there is not -if, there is not -of
      there_is_field = .false.
    elseif (.not. is_arg('-of')) then
      !there is -if, there is not -of
      call set( infield, get_post_arg('-if'), 1, fit=.true.)
      p = index(infield(1), '.', back=.true.)
      call set(outfield, trim(outmesh)//'__'//trim(infield(1)(1:p-1))//'.'//trim(FLDB(id_mesh_ext(outext))%field_ext), 1, &
      &fit=.true.)
    else
      !there is not -if, there is -of
      call error('(module_feconv/fe_conv) option -fi is mandatory to read external fields')
    end if
  else
    !infield is mesh external, outfield is mesh internal
    if (is_arg('-if') .and. is_arg('-of')) then
      !there is -if, there is -of
      call set( infield, get_post_arg('-if'), 1, fit=.true.)
      call set(outfield, get_post_arg('-of'), 1, fit=.true.)
    elseif (.not. is_arg('-if') .and. .not. is_arg('-of')) then
      !there is not -if, there is not -of
      there_is_field = .false.
    elseif (.not. is_arg('-of')) then
      !there is -if, there is not -of
      call set( infield, get_post_arg('-if'), 1, fit=.true.)
      p = index(infield(1), '.', back=.true.); call set(outfield, trim(infield(1)(1:p-1)), 1, fit=.true.)
    else
      !there is not -if, there is -of
      call error('(module_feconv/fe_conv) option -if is mandatory to read external fields')
    end if
  end if
elseif (FLDB(id_mesh_ext(outext))%is_field_outside) then
  !infield is mesh internal, outfield is mesh external
  if (is_arg('-if') .and. is_arg('-of')) then
    !there is -if, there is -of
    call set( infield, get_post_arg('-if'), 1, fit=.true.)
    call set(outfield, get_post_arg('-of'), 1, fit=.true.)
  elseif (.not. is_arg('-if') .and. .not. is_arg('-of')) then
    !there is not -if, there is not -of
    there_is_field = .false.
  elseif (.not. is_arg('-of')) then
    !there is -if, there is not -of
    call set( infield, get_post_arg('-if'), 1, fit=.true.)
    call set(outfield, trim(outmesh)//'__'//trim(infield(1))//'.'//trim(FLDB(id_mesh_ext(outext))%field_ext), 1, fit=.true.)
  else
    !there is not -if, there is -of
    call set(infield, '*', 1, fit=.true.)
    call set(outfield, get_post_arg('-of'), 1, fit=.true.)
  end if
else
  !infield and outfield are both mesh internal
  if (is_arg('-if') .and. is_arg('-of')) then
    !there is -if, there is -of
    call set( infield, get_post_arg('-if'), 1, fit=.true.)
    call set(outfield, get_post_arg('-of'), 1, fit=.true.)
  elseif (.not. is_arg('-if') .and. .not. is_arg('-of')) then
    !there is not -if, there is not -of
    call set( infield, '*', 1, fit=.true.)
    call set(outfield, '*', 1, fit=.true.)
  elseif (.not. is_arg('-of')) then
    !there is -if, there is not -of
    call set( infield, get_post_arg('-if'), 1, fit=.true.)
    call set(outfield, trim(infield(1)), 1, fit=.true.)
  else
    !there is not -if, there is -of
    call error('(module_feconv/fe_conv) option -fi is mandatory to read a specific field.')
  end if
end if
  
!si es oM, load aparte; si es iM, se engorda load
!escritura, lo mismo

!read mesh
is_pmh = .false.
select case (trim(lcase(adjustlt(inext))))
case('mfm')
  print '(a)', 'Loading MFM mesh file...'
  call load_mfm(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('mum')
  print '(a)', 'Loading MUM mesh file...'
  call load_mum(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('msh')
  print '(a)', 'Loading ANSYS mesh file...'
  call load_msh(infile, pmh); is_pmh = .true.
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
  call load_mphtxt(infile, pmh); is_pmh = .true.
  print '(a)', 'Done!'
case('pf3')
  print '(a)', 'Loading FLUX mesh file...'
  call load_pf3(infile, pmh); is_pmh = .true.
  print '(a)', 'Done!'
case('vtu')
  print '(a)', 'Loading MFM mesh file...'
  call load_vtu( infile, pmh); is_pmh = .true.
  print '(a)', 'Done!'
case default
  call error('(module_feconv/fe_conv) input file extension not implemented: '//trim(adjustlt(inext)))
end select

!transform
if (is_arg('-l1')) then
  print '(/a)', 'Converting mesh into Lagrange P1 mesh...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call to_l1(pmh); is_pmh = .true.
elseif (is_arg('-l2')) then
  print '(/a)', 'Converting Lagrange P1 mesh into Lagrange P2 mesh...'
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call lagr2l2(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm); is_pmh = .false.
elseif (is_arg('-rt')) then
  print '(/a)', 'Converting Lagrange P1 mesh into Raviart-Thomas (face) mesh...'
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call lagr2rt(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm); is_pmh = .false.
elseif (is_arg('-nd')) then
  print '(/a)', 'Converting Lagrange P1 mesh into Whitney (edge) mesh...'
  if(is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call lagr2nd(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm); is_pmh = .false.
end if

!bandwidth optimization
if (is_arg('-cm')) then
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call cuthill_mckee(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, z); is_pmh = .false.
end if

!save mesh
select case (trim(adjustlt(outext)))
case('mfm')
  print '(/a)', 'Saving MFM mesh file...'
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call save_mfm(outfile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
case('mum')
  print '(/a)', 'Saving MUM mesh file...'
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call save_mum(outfile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  print '(a)', 'Done!'
!case('vtu')
!  print '(/a)', 'Saving VTU mesh file...'
!  if(allocated(pmh%pc)) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
!  call save_vtu(outfile, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
!  print '(a)', 'Done!'
case('vtu')
  print '(/a)', 'Saving VTU mesh file...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_vtu2(outfile, pmh)
  print '(a)', 'Done!'
case('mphtxt')
  print '(/a)', 'Saving COMSOL mesh file...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_mphtxt(outfile, pmh)
  print '(a)', 'Done!'
case('unv')
  print '(/a)', 'Saving I-DEAS UNV mesh file...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_unv(outfile, get_unit(), pmh)
  print '(a)', 'Done!'
case('pf3')
  print '(/a)', 'Saving FLUX mesh file...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_pf3(outfile, pmh)
  print '(a)', 'Done!'
case('msh')
  print '(/a)', 'Saving ANSYS mesh file...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_msh(outfile, pmh)
  print '(a)', 'Done!'
case('mesh')
  print '(/a)', 'Saving FreFem++ mesh file...'
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  call save_freefem(outfile, get_unit(), pmh)
  print '(a)', 'Done!'
end select !case default, already checked before reading infile
end subroutine

end module
