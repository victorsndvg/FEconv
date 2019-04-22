module module_feconv
!-----------------------------------------------------------------------
! Module to convert between several mesh and FE field formats
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
!
! PUBLIC PROCEDURES:
! convert: converts between several mesh and FE field formats
! is_arg: returns true when the argument is present
!-----------------------------------------------------------------------
use basicmod, only: real64, maxpath, slash, error, adjustlt, lcase, word_count, get_unit, set, get_arg, is_arg, &
                    get_post_arg, args_count, set_args, file_exists, operator(.IsNewerThan.), report_option
use module_transform_fcnv, only: lagr2l2, lagr2rt, lagr2nd, lagr2nd2, to_l1
use module_cuthill_mckee_fcnv, only: cuthill_mckee
use module_msh_fcnv, only: load_msh,save_msh
use module_unv_fcnv, only: load_unv,save_unv
use module_patran_fcnv, only: load_patran
use module_mfm_fcnv, only: load_mfm, save_mfm
use module_mum_fcnv, only: load_mum, save_mum
use module_vtu_fcnv, only: load_vtu, save_vtu, type_cell
use module_pvd_fcnv, only: load_pvd, save_pvd
use module_mphtxt_fcnv, only: load_mphtxt,save_mphtxt
use module_pf3_fcnv, only: load_pf3,save_pf3
!use module_tra, only: load_tra,save_tra
use module_field_database_fcnv, only: FLDB, id_mesh_ext, id_field_ext
use module_mff_fcnv, only: load_mff, save_mff
use module_muf_fcnv, only: load_muf, save_muf
use module_freefem_fcnv, only: save_freefem_msh, save_freefem_mesh, load_freefem_msh, load_freefem_mesh
use module_pmh_fcnv
use module_fem_extract_fcnv, only: extract_mesh, extract_ref
use module_gmsh_fcnv, only: load_gmsh, save_gmsh
use module_dex_fcnv, only: load_dex, save_dex
use module_ip_fcnv, only: load_ip, save_ip
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

!Private procedures
private :: get_fieldfile, get_fieldname

contains

!-----------------------------------------------------------------------
! convert: converts between several mesh and FE field formats
!-----------------------------------------------------------------------
subroutine convert(argstr, inpmh, outpmh, force)
character(*),   optional, intent(in)    :: argstr !! String with arguments.
type(pmh_mesh), optional, intent(inout) :: inpmh  !! Input PMH structure.
type(pmh_mesh), optional, intent(inout) :: outpmh !! Output PMH structure.
logical,        optional, intent(in)    :: force  !! Whether to force to save `outfile`.
character(maxpath) :: infile=' ', inmesh=' ', inext=' ', outfile=' ', outmesh=' ', outext=' '
character(maxpath) :: infext=' ', outfext=' ', outpath = ' '!,fieldfilename = ' '
character(maxpath), allocatable :: infieldfile(:), outfieldfile(:),infieldname(:), outfieldname(:)
character(maxpath) :: str
integer :: p, nargs, q, comp
integer, allocatable :: nsd0(:), chref(:)
logical :: there_is_field, force_to_save
!Variables for extratction
integer,      allocatable :: submm(:,:), subnrv(:,:), subnra(:,:), subnrc(:,:), subnsd(:), globv(:), globel(:)
real(real64), allocatable :: subz(:,:)
real(real64)              :: padval

!call report_option('info', 'std')
if(present(argstr)) then
  call set_args(argstr)
else
  call info('String options for convert not found; reading options from command line.')
endif

!find infile and outfile at the end of the arguments
nargs = args_count()

if(present(inpmh)) then
  ! with inpmh present, only outfile is read as the last argument
  outfile = get_arg(nargs);     p      = index( outfile, '.', back=.true.)
  outmesh = outfile(1:p-1);     outext = lcase(outfile(p+1:len_trim(outfile)))
elseif(present(outpmh) .or. is_arg('-l')) then
  ! with outpmh present or with argument -l, only infile is read as the last argument
  infile = get_arg(nargs);      p = index(infile, '.', back=.true.)
  inmesh = infile(1:p-1);       inext = lcase(infile(p+1:len_trim(infile)))
else
  ! otherwise, inmesh outmesh are read as the two last arguments
  infile = get_arg(nargs-1);    p      = index( infile, '.', back=.true.)
  inmesh = infile(1:p-1);       inext  =  lcase(infile(p+1:len_trim( infile)))
  outfile = get_arg(nargs);     p      = index(outfile, '.', back=.true.)
  outmesh = outfile(1:p-1);     outext = lcase(outfile(p+1:len_trim(outfile)))
  p = index(outfile, slash(), back=.true.); outpath = outfile(1:p)
endif

!check mesh names and extensions
if (.not. present(inpmh)) then
  if (len_trim(infile)  == 0)    call error('(module_feconv::convert) unable to find input file: '//trim(infile))
  if (.not. file_exists(infile)) call error('(module_feconv::convert) input file does not exist: '//trim(infile))
  if (len_trim(inext)   == 0)    call error('(module_feconv::convert) unable to find input file extension: '//trim(inext))
end if
if(.not. (present(outpmh) .or. is_arg('-l'))) then
  if (len_trim(outfile) == 0)    call error('(module_feconv::convert) unable to find output file: '//trim(outfile))
  if (len_trim(outext)  == 0)    call error('(module_feconv::convert) unable to find output file extension: '//trim(outext))
  select case (trim(adjustlt(outext))) !check outfile extension now (avoid reading infile when outfile is invalid)
  case('mfm', 'mum', 'vtu', 'mphtxt', 'unv', 'pf3', 'msh', 'mesh', 'pmh', 'pvd')
    continue
  case default
    call error('(module_feconv::convert) output file extension not implemented: '//trim(adjustlt(outext)))
  end select
endif

! set force_to_save
force_to_save = .true.
if (present(force)) force_to_save = force
if (.not. present(outpmh) .and. .not. is_arg('-l')) then 
  if (.not. file_exists(outfile)) force_to_save = .true.
end if
if (.not. present(inpmh) .and. .not. present(outpmh) .and. .not. is_arg('-l')) then 
  if (infile .IsNewerThan. outfile) force_to_save = .true.
end if
call info('Forcing to save output file: '//string(force_to_save))

!check isoparametric option, for UNV only
!if (trim(adjustlt(inext)) /= 'unv' .and. is_arg('-is')) call error('(module_feconv/fe_conv) only UNV input files can '//&
!&'manage -is option.')
!options for mesh transformation (-l1, -l2, -rt, -nd, -nd2 and -cm) are incompatible with fields (-if, -of)
if ( (is_arg('-l1') .or. is_arg('-l2') .or. is_arg('-rt') .or. is_arg('-nd') .or. is_arg('-nd2') .or. is_arg('-cm')) .and. &
     (is_arg('-if') .or. is_arg('-of')) ) call error('(module_feconv/fe_conv) options for mesh transformation (-l1, -l2, '//&
     &'-rt, -nd, -nd2 and -cm) are incompatible with fields (-if, -of).')
!set PMH mesh tolerance (all load procedures must consider intent(inout) for PMH argument)
if (is_arg('-t')) then
  pmh%ztol = dble(get_post_arg('-t'))
end if

! field selection
there_is_field = .true.
if (present(outpmh)) then
  ! inext does not exist but -if must be read
  if (is_arg('-if')) then
    call get_fieldfile('-if', infieldfile, infext)
    if (is_arg('-in')) call get_fieldname('-in', infieldname)
  else
    there_is_field = .false.
  end if
elseif (is_arg('-l')) then
  ! outext does not exist but -if must be read
  if (is_arg('-if')) then
    call get_fieldfile('-if', infieldfile, infext)
    if (is_arg('-in')) call get_fieldname('-in', infieldname)
  else
    there_is_field = .false.
  endif
elseif (present(inpmh)) then
  ! outext does not exist but -of must be read
  if (is_arg('-of')) then
    call get_fieldfile('-of', outfieldfile, outfext)
    if (is_arg('-on')) call get_fieldname('-on', outfieldname)
  else
    there_is_field = .false.
  endif
else
  ! both inext and outext exist
  if (FLDB(id_mesh_ext(inext))%is_field_outside) then
    if (FLDB(id_mesh_ext(outext))%is_field_outside) then 
      ! infieldfile and outfieldfile are both mesh external
      if (.not. is_arg('-if') .and. is_arg('-of')) then
        call error('(module_feconv/fe_conv) option -if is mandatory to read external fields')
      elseif (.not. is_arg('-if')) then                  ! there is neither -if nor -of
        there_is_field = .false.
      elseif (is_arg('-if') .or. is_arg('-of')) then     ! there is -if or -of
        if (is_arg('-if')) then
          call get_fieldfile('-if', infieldfile, infext)
          ! read -in when field extension is 'ip'
          if (id_field_ext(infext) == id_field_ext('ip') .and. is_arg('-in')) call get_fieldname('-in', infieldname)
        endif
        if (is_arg('-of')) then
          call get_fieldfile('-of', outfieldfile, outfext)
          ! read -on when field extension are 'ip' or 'dex'
          if (is_arg('-on') .and. (id_field_ext(outfext) == id_field_ext('ip') .or. &
                                   id_field_ext(outfext) == id_field_ext('dex'))) call get_fieldname('-on', outfieldname)
        endif
      end if
    else
      ! infieldfile is mesh external, outfieldfile is mesh internal
      if (is_arg('-if')) then
        call get_fieldfile('-if', infieldfile, infext)
        ! read -in when field extension 'ip' or 'mff'
        if (is_arg('-in') .and. (id_field_ext(infext) == id_field_ext('ip') .or. &
                                 id_field_ext(infext) == id_field_ext('mff'))) call get_fieldname('-in', infieldname)
        if (is_arg('-on')) call get_fieldname('-on', outfieldname)
      else
        there_is_field = .false.
      end if
    end if
  elseif (FLDB(id_mesh_ext(outext))%is_field_outside) then
    ! infieldfile is mesh internal, outfieldfile is mesh external
    if (is_arg('-of')) then
      call get_fieldfile('-of', outfieldfile, outfext)
      ! read -in when mesh extension is not 'pf3'
      if (is_arg('-in') .and. id_mesh_ext(inext) /= id_mesh_ext('pf3')) call get_fieldname('-in', infieldname) 
      ! read -on when field extension is neither 'mff' nor 'muf'
       if (is_arg('-on') .and. (id_field_ext(inext) /= id_field_ext('mff') .and. &
                                id_field_ext(inext) /= id_field_ext('muf'))) call get_fieldname('-on', outfieldname)
    end if
  else
    ! infieldfile and outfieldfile are both mesh internal
    if (is_arg('-in')) call get_fieldname('-in', infieldname)
    if (is_arg('-on')) call get_fieldname('-on', outfieldname)
  end if
endif

! Sets the field padding value
if (is_arg('-pad')) then
  padval = dble(get_post_arg('-pad'))
else
  padval = 0._real64
endif

!read mesh
if (present(inpmh)) then
  pmh = inpmh
  is_pmh = .true.
else
  is_pmh = .false.
  select case (trim(lcase(adjustlt(inext))))
  case('mfm')
    call info('Loading MFM mesh file...')
    call load_mfm(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
      call info('Done!')
  case('mum')
    call info('Loading MUM mesh file...')
    call load_mum(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    call info('Done!')
  case('msh')
    if (is_arg('-ff')) then !FreeFem++
      call info('Loading FreFem++ (.msh) mesh file...')
      call load_freefem_msh(infile, get_unit(), pmh); is_pmh = .true.
    elseif (is_arg('-gm')) then !Gmsh
      call info('Loading Gmsh (.msh) mesh file...')
      call load_gmsh(infile, get_unit(), pmh); is_pmh = .true.
    else !ANSYS
      call info('Loading ANSYS mesh file...')
      call load_msh(infile, pmh); is_pmh = .true.
    end if
    call info('Done!')
  case('unv')
    call info('Loading UNV mesh file...')
    call load_unv(infile, pmh, padval, infieldname, is_arg('-ca')); is_pmh = .true.
    call info('Done!')
  case('bdf')
    call info('Loading MD Nastran input file...')
    call load_patran(infile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    call info('Done!')
  case('mphtxt')
    call info('Loading COMSOL mesh file...')
    call load_mphtxt(infile, pmh); is_pmh = .true.
    call info('Done!')
  case('pf3')
    call info('Loading FLUX mesh file...')
    call load_pf3(infile, pmh); is_pmh = .true.
    call info('Done!')
  case('vtu')
    call info('Loading VTU mesh file...')
    call load_vtu( infile, pmh, infieldname); is_pmh = .true.
    call info('Done!')
  case('pvd')
    call info('Loading PVD file...')
    call load_pvd( infile, pmh, infieldname); is_pmh = .true.
    call info('Done!')
  case('mesh')
    call info('Loading FreFem++ (Tetrahedral Lagrange P1) MESH file...')
    call load_freefem_mesh(infile, get_unit(), pmh); is_pmh = .true.
    call info('Done!')
  case default
    call error('(module_feconv/fe_conv) input file extension not implemented: '//trim(adjustlt(inext)))
  end select
end if

! Read field files
if (there_is_field .and. is_arg('-if') .and. (.not. present(inpmh))) then
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  is_pmh = .true.
  select case (trim(lcase(adjustlt(infext))))
    case('mff')
      call load_mff(pmh, infieldfile, infieldname)
    case('muf')
      call load_muf(pmh, infieldfile, infieldname)
    case('dex')
      call load_dex(pmh, infieldfile, infieldname)
    case('ip')
      call load_ip(pmh, infieldfile, infieldname, outfieldname)
  end select
endif

! Show PMH info in screen
if(is_arg('-l')) then
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  is_pmh = .true.
  call save_pmh(pmh, ' ', with_values=.false.) !only printed in screen, filename is not required
  stop
endif

! Remove a component of the space dimension
if (is_arg('-rc')) then
 comp = int(get_post_arg('-rc'))
 if (.not. is_pmh) then
   call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
   is_pmh=.true.
 endif
 call remove_coordinate(pmh, comp)
endif

!extract (only for Lagrange P1 meshes)
if (is_arg('-es')) then
  str = get_post_arg('-es')
  call info('Extracting subdomain(s) '//trim(str)//'...')
  p = index(str, '[')
  if (p == 0) then !a single subdomain ref.
    call set(nsd0, int(str), 1, fit=.true.)
  else !several subdomain refs. enclosed in [] and separated by ,
    q = index(str, ']', back=.true.)
    call alloc(nsd0, word_count(str(p+1:q-1),','))
    read(str(p+1:q-1),*) nsd0
  end if
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd); is_pmh = .false.
  if (nver /= nnod) call error('(module_feconv/fe_conv) extraction is only available for Lagrange P1 meshes.')
  call extract_mesh(nver, mm, z, nsd, nsd0, submm, subz, globv, globel)
  call extract_ref(nrv, nra, nrc, nsd, subnrv, subnra, subnrc, subnsd, globel)
  nel  = size(submm, 2)
  nver = size(subz,  2)
  nnod = nver
  call move_alloc(from=submm,  to=mm)
  call move_alloc(from=subnrv, to=nrv)
  call move_alloc(from=subnra, to=nra)
  call move_alloc(from=subnrc, to=nrc)
  call move_alloc(from=subz,   to=z)
  call move_alloc(from=subnsd, to=nsd)
  call info('Done!')
end if

!transform
if (is_arg('-l1')) then
  call info('Converting mesh into Lagrange P1 mesh...')
  if (.not. is_pmh) then
    call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    is_pmh = .true.
  end if
  call to_l1(pmh)
  call info('Done!')
elseif (is_arg('-l2')) then
  call info('Converting Lagrange P1 mesh into Lagrange P2 mesh...')
  if (is_pmh) then
    call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    is_pmh = .false.
  end if
  call lagr2l2(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
  call info('Done!')
elseif (is_arg('-rt')) then
  call info('Converting Lagrange mesh into Raviart-Thomas (face) mesh...')
  if (is_pmh) then
    call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    is_pmh = .false.
  end if
  call lagr2rt(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
  call info('Done!')
elseif (is_arg('-nd')) then
  call info('Converting Lagrange mesh into Whitney (edge) mesh...')
  if (is_pmh) then
    call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    is_pmh = .false.
  end if
  call lagr2nd(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
  call info('Done!')
elseif (is_arg('-nd2')) then
  call info('Converting Lagrange mesh into Whitney order 2 (edge) mesh...')
  if (is_pmh) then
    call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    is_pmh = .false.
  end if
  call lagr2nd2(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
  call info('Done!')
end if

!bandwidth optimization
if (is_arg('-cm')) then
  if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
  call cuthill_mckee(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, z); is_pmh = .false.
end if

!cell to node
if (is_arg('-cn')) then
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  is_pmh = .true.
  call cell2node(pmh)
end if

!change references
if (is_arg('-ch')) then
  str = get_post_arg('-ch')
  p = index(str, '[')
  if (p > 0) then
    q = index(str, ']', back=.true.)
    call alloc(chref, word_count(str(p+1:q-1),','))
    read(str(p+1:q-1),*) chref
    if (mod(size(chref,1), 2) == 0) call error('(module_feconv/convert) argument after -ch must be: [<topological-value>, '//&
    &'old-ref_1, new-ref_1, ..., old-ref_n, new-ref_n].')
  else
    call error('(module_feconv/convert) argument after -ch must be enclosed in square brackets.')
  end if
  if (.not.is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  is_pmh = .true.
  call change_pmh_references(pmh, chref)
end if

!save mesh
if(present(outpmh)) then
  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  outpmh = pmh
elseif (force_to_save) then
  select case (trim(adjustlt(outext)))
  case('mfm')
    call info('Saving MFM mesh file...')
    if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    call save_mfm(outfile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    call info('Done!')
  case('mum')
    call info('Saving MUM mesh file...')
    if (is_pmh) call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    call save_mum(outfile, get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    call info('Done!')
  case('vtu')
    call info('Saving VTU mesh file...')
    if (is_pmh) then
      call save_vtu(outfile, pmh, infieldname, outfieldname, padval)
    else
      call save_vtu(outfile, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
    endif
    call info('Done!')
  !case('vtu')
  !  print '(/a)', 'Saving VTU mesh file...'
  !  if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
  !  call save_vtu2(outfile, pmh)
  !  call info('Done!')
  case('pvd')
    call info('Saving PVD file...')
    if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    call save_pvd(outfile, pmh,infieldname, outfieldname, padval)
    call info('Done!')
  case('mphtxt')
    call info('Saving COMSOL mesh file...')
    if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    call save_mphtxt(outfile, pmh)
    call info('Done!')
  case('unv')
    call info('Saving I-DEAS UNV mesh file...')
    if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    call save_unv(outfile, get_unit(), pmh, infieldname, outfieldname, is_arg('-ca'))
    call info('Done!')
  case('pf3')
    call info('Saving FLUX mesh file...')
    if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    call save_pf3(outfile, pmh, infieldname, outfieldname, outpath)
    call info('Done!')
  case('msh')
    if (is_arg('-ff')) then !FreeFem++
      call info('Saving FreFem++ mesh file...')
      if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
      call save_freefem_msh(outfile, get_unit(), pmh)
    elseif (is_arg('-gm')) then !Gmsh
      call info('Saving Gmsh mesh file...')
      if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
      call save_gmsh(outfile, get_unit(), pmh)
    else !ANSYS
      call info('Saving ANSYS mesh file...')
      if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
      call save_msh(outfile, pmh)
    end if
    call info('Done!')
  case('mesh')
    call info('Saving FreFem++ mesh file...')
    if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    call save_freefem_mesh(outfile, get_unit(), pmh)
    call info('Done!')
  case('pmh')
    call info('Saving PMH mesh file...')
    if (.not. is_pmh) call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
    call save_pmh(pmh, outfile)
    call info('Done!')
  end select !case default, already checked before reading infile
  !save fields
  if(is_pmh .and. there_is_field .and. is_arg('-of').and. (.not. present(outpmh))) then
    select case (trim(lcase(adjustlt(outfext))))
      case('mff')
        call save_mff(pmh, outfieldfile, outpath)
      case('muf')
        call save_muf(pmh, outfieldfile, outpath)
      case('dex')
        call save_dex(pmh, infieldname, outfieldname, outfieldfile)
      case('ip')
        call save_ip(pmh, outfieldfile, infieldname, outfieldname)
      case default
        call info('Field file extension "'//trim(lcase(adjustlt(outfext)))//'" not supported!')
    end select
  endif
endif
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! get_fieldfile: get variables related with option -if and -of
!-----------------------------------------------------------------------
subroutine get_fieldfile(opt, fieldfile, fext)
character(*), intent(in)                 :: opt
character(*), intent(inout), allocatable :: fieldfile(:)
character(*), intent(inout)              :: fext
character(maxpath) :: str
integer :: p, q

str = get_post_arg(trim(adjustl(opt)))
p = index(str, '[')
if (p == 0) then !a single field name
  call set(fieldfile, str, 1, fit=.true.)
else !several subdomain refs. enclosed in [] and separated by ,
  q = index(str, ']', back=.true.)
  call alloc(fieldfile, word_count(str(p+1:q-1),','))
  read(str(p+1:q-1),*) fieldfile
end if
p = index(fieldfile(1), '.', back=.true.)
fext =  fieldfile(1)(p+1:len_trim(fieldfile(1)))
end subroutine

!-----------------------------------------------------------------------
! get_fieldname: get variables related with option -in and -on
!-----------------------------------------------------------------------
subroutine get_fieldname(opt, fieldname)
character(*), intent(in)                 :: opt
character(*), intent(inout), allocatable :: fieldname(:)
character(maxpath) :: str
integer :: p, q

str = get_post_arg(trim(adjustl(opt)))
p = index(str, '[')
if (p == 0) then !a single field name
  call set(fieldname, str, 1, fit=.true.)
else !several subdomain refs. enclosed in [] and separated by ,
  q = index(str, ']', back=.true.)
  call alloc(fieldname, word_count(str(p+1:q-1),','))
  read(str(p+1:q-1),*) fieldname
end if
end subroutine

end module
