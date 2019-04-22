module module_pmh_fcnv
!-----------------------------------------------------------------------
! Module to manage piecewise meshes (PMH)
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 19/01/2014
!
! PUBLIC PROCEDURES:
!   save_pmh: save mesh into a PMH file
!   pmh2mfm: convert a PMH structure into a MFM one
!   mfm2pmh: convert a MFM mesh into a PMH structure
!   build_vertices: build vertex connectivities and coordinates from node information (P1, P2 only)
!   build_node_coordinates: build node coordinates from vertex information
!   cell2node: calculate a node field form a cell one
!   get_field_num_shots: returns the number of shots giving a field name
!   get_num_shots: returns an array with the number of shots of all fields
!   get_piece_max_top_dim: returns the maximum topological dimension
!   remove_coordinate: reduces the space dimension of the mesh removing the chosen coordinate
!   change_pmh_references: changes pmh references
!
! REMARKS:
!   A mesh is divided into pieces
!   Each piece has common vertices/nodes and it can contain one or several groups of elements
!   Each group of elements belong to a type of element defined in module_fe_database_pmh
!   In the group(s) with maximal topological dimension:
!     - Vertex connectivities (mm) are mandatory;
!     - Node conectivities (nn) are only required for non Lagrange-P1 elements
!   In the group(s) with topological dimension smaller than the maximal, vertex connectivities (mm) are relevant
!   Variable z always contains vertex coordinates (nodes coordinates can be deduced from z and the element type)
!   Global numbering starts at 1 for elements, vertices and nodes
!-----------------------------------------------------------------------
use basicmod, only: real64, output_unit, maxpath, error, info, string, int, dble, alloc, dealloc, set, reduce, sfind, &
                    find_first, find_sorted, sort, insert, insert_sorted, bsearch, is_arg, get_post_arg, feed, empty, &
                     unique, get_unit, sunique
use module_fe_database_pmh_fcnv, only : FEDB, check_fe, VF_WEDG
implicit none

!Types

type field
  character(maxpath)        :: name
  real(real64), allocatable :: param(:)   !nshot
  real(real64), allocatable :: step(:)    !nshot
  real(real64), allocatable :: val(:,:,:) !ncomp x nnod x nshot
end type

type elgroup
  integer                  :: nel  = 0 !total number of elements
  integer                  :: type = 0 !element type (one of those defined in module_eltype)
  integer,     allocatable :: nn(:,:)  !global numbering of nodes
  integer,     allocatable :: mm(:,:)  !global numbering of vertices
  integer,     allocatable :: ref(:)   !reference numbering
  type(field), allocatable :: fi(:)    !fields on elements
end type

type piece
  integer                    :: nnod = 0 !total number of nodes
  integer                    :: nver = 0 !total number of vertices
  integer                    :: dim  = 0 !space dimension of the node/vertex coordinates
  real(real64),  allocatable :: z(:,:)   !vertex coordinates
  type(elgroup), allocatable :: el(:)    !element groups
  type(field),   allocatable :: fi(:)    !fields on nodes
end type

type pmh_mesh
  type(piece), allocatable :: pc(:) !pieces that compose the mesh
  real(real64)             :: ztol = epsilon(0._real64) !mesh tolerance
end type

!Constants
!Tetrahedra used to calculate a sufficient condition to ensure positive Jacobian for an hexahedron
!(see http://www.math.udel.edu/~szhang/research/p/subtettest.pdf)
integer, parameter, private :: Pc(32,4) = reshape([ &
1, 2, 3, 7,  2, 3, 4, 8,  3, 4, 1, 5,  4, 1, 2, 6,  5, 8, 7, 3,  8, 7, 6, 2,  7, 6, 5, 1,  6, 5, 8, 4,  &
2, 1, 5, 8,  4, 3, 7, 6,  1, 4, 8, 7,  3, 2, 6, 5,  1, 2, 3, 5,  2, 3, 4, 6,  3, 4, 1, 7,  4, 1, 2, 8,  &
5, 8, 7, 1,  8, 7, 6, 4,  7, 6, 5, 3,  6, 5, 8, 2,  5, 1, 4, 6,  3, 7, 8, 2,  4, 8, 5, 3,  2, 6, 7, 1,  &
1, 2, 4, 5,  2, 3, 1, 6,  3, 4, 2, 7,  4, 1, 3, 8,  5, 8, 6, 1,  6, 5, 7, 2,  7, 6, 8, 3,  8, 7, 5, 4], [32, 4], order=[2,1])

!Private procedures
private :: swap, reorder_nodes_element_P2, QJ_pos, detDFT_pos, reorder_nodes_P2, positive_jacobian, &
           cell2node_real_pmh

interface cell2node;     module procedure cell2node_real_pmh;     end interface

contains

!-----------------------------------------------------------------------
! save_pmh: save pmh
!-----------------------------------------------------------------------
!subroutine save_pmh(filename, pmh, with_values)
!character(*),   intent(in) :: filename !mesh filename
!type(pmh_mesh), intent(in) :: pmh      !pmh structure
!logical, optional          :: with_values
!integer                    :: iu       !file unit
!logical :: wv
!integer :: i, j, k, ip, ig, ios
! 
!wv = .true.
!if (present(with_values)) wv = with_values
! 
!if(wv) then
!  iu = get_unit()
!  open (unit=iu, file=filename, form='formatted', position='rewind', iostat=ios)
!  if (ios /= 0) call error('module_pmh/save_pmh: open error #'//trim(string(ios)))
!else
!  iu = output_unit
!endif
!write(iu, '(a/)') '<?xml version="1.0" encoding="UTF-8" ?>'
!write(iu, '(a)') '<pmh>'
!do ip = 1, size(pmh%pc,1)
!  write(iu, '(2x,a)') '<piece name="'//trim(string(ip))//'">'
!  write(iu, '(4x,a)') '<nnod> '//trim(string(pmh%pc(ip)%nnod))//' </nnod>'
!  write(iu, '(4x,a)') '<nver> '//trim(string(pmh%pc(ip)%nver))//' </nver>'
!  write(iu, '(4x,a)')  '<dim> '//trim(string(pmh%pc(ip)%dim))//' </dim>'
!  if(wv) then
!    write(iu, '(4x,a)') '<z>'
!    do k = 1, pmh%pc(ip)%nver; do j = 1, pmh%pc(ip)%dim; call feed(iu, string(pmh%pc(ip)%z(j,k))); end do; end do; call empty(iu)
!    write(iu, '(4x,a)') '</z>'
!  endif
!  do ig = 1, size(pmh%pc(ip)%el,1)
!    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
!      write(iu, '(4x,a)') '<element_group name="'//trim(string(ig))//'">'
!      write(iu, '(6x,a)')  '<nel> '//trim(string(elg%nel))//' </nel>'
!      write(iu, '(6x,a)') '<type> '//trim(string(elg%type))//' </type>'
!      write(iu, '(6x,a)') '<desc>'
!      write(iu,    '(a)') trim(FEDB(elg%type)%desc)
!      write(iu, '(6x,a)') '</desc>'
! 
!      if(wv) then
!        if (allocated(elg%nn)) then
!          write(iu, '(6x,a)') '<nn>'
!          do k = 1, elg%nel; do i = 1, FEDB(tp)%lnn; call feed(iu, string(elg%nn(i,k))); end do; end do; call empty(iu)
!          write(iu, '(6x,a)') '</nn>'
!        end if
!        write(iu, '(6x,a)') '<mm>'
!        do k = 1, elg%nel; do i = 1, FEDB(tp)%lnv; call feed(iu, string(elg%mm(i,k))); end do; end do; call empty(iu)
!        write(iu, '(6x,a)') '</mm>'
!        write(iu, '(6x,a)') '<ref>'
!        do k = 1, elg%nel; call feed(iu, string(elg%ref(k))); end do; call empty(iu)
!        write(iu, '(6x,a)') '</ref>'
!      endif
!      write(iu, '(4x,a)') '</element_group>'
!    end associate
!  end do
!  write(iu, '(2x,a)') '</piece>'
!end do
!write(iu, '(a)') '</pmh>'
!if(wv) close(iu)
!end subroutine

!-----------------------------------------------------------------------
! save_pmh: save pmh
!-----------------------------------------------------------------------
subroutine save_pmh(pmh, filename, with_values)
!! Prints or saves a PMH structure.  
!! 
!! When argument `with_values' is true or not present, the PMH structure is completely saved in `filename`; 
!! otherwise, the scalar data of the structure is printed (and `filename` is not used).  
!!
!! @note The renovated version of this procedure takes into account the PMH fields. It also includes scalar data as attibutes in 
!! tags  `piece`, `field`, `element_group`. For example,  
!! `<piece name="1" nnod="2236" nver="2155" dim="3">`  
type(pmh_mesh),         intent(in) :: pmh         !! Input PMH structure
character(*), optional, intent(in) :: filename    !! Output filename.
logical,      optional, intent(in) :: with_values !! Whether the arrays values are preinted or not.
integer :: iu
logical :: wv
integer :: i, j, k, ip, ig, ifi, ios

wv = .true.
if (present(with_values)) wv = with_values

iu = output_unit
if (wv) then
  iu = get_unit()
  if (.not. present(filename)) call error('(module_pmh::save_pmh) filename  not present.')
  open (unit=iu, file=filename, form='formatted', position='rewind', iostat=ios)
  if (ios /= 0) call error('module_pmh/save_pmh: open error #'//trim(string(ios)))
endif
write(iu, '(a/)') '<?xml version="1.0" encoding="UTF-8" ?>'
write(iu, '(a)') '<pmh>'
do ip = 1, size(pmh%pc,1)
  write(iu, '(2x,a)') '<piece name="'//trim(string(ip))//&
                         & '" nnod="'//trim(string(pmh%pc(ip)%nnod))//&
                         & '" nver="'//trim(string(pmh%pc(ip)%nver))//&
                         & '" dim="'//trim(string(pmh%pc(ip)%dim))//'">'
  if(wv) then
    write(iu, '(4x,a)') '<z>'
    do k = 1, pmh%pc(ip)%nver; do j = 1, pmh%pc(ip)%dim; call feed(iu, string(pmh%pc(ip)%z(j,k))); end do; end do; call empty(iu)
    write(iu, '(4x,a)') '</z>'
  endif
  if(allocated(pmh%pc(ip)%fi)) then
    do ifi=1, size(pmh%pc(ip)%fi,1)
      if(.not. allocated(pmh%pc(ip)%fi(ifi)%val)) cycle
      write(iu, '(4x,a)') '<field name="'//trim(pmh%pc(ip)%fi(ifi)%name)//&
                             & '" ncomp="'//trim(string(size(pmh%pc(ip)%fi(ifi)%val,1)))//&
                             & '" nshot="'//trim(string(size(pmh%pc(ip)%fi(ifi)%param,1)))//'">'
      write(iu, '(4x,a)') '</field>'
    enddo
  endif
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
      write(iu, '(4x,a)') '<element_group name="'//trim(string(ig))//&
                                     & '" nel="'//trim(string(elg%nel))//&
                                     & '" type="'//trim(string(elg%type))//'">'
      write(iu, '(6x,a)') '<desc> '//trim(FEDB(elg%type)%desc)//' </desc>'
      if(wv) then
        if (allocated(elg%nn)) then
          write(iu, '(6x,a)') '<nn>'
          do k = 1, elg%nel; do i = 1, FEDB(tp)%lnn; call feed(iu, string(elg%nn(i,k))); end do; end do; call empty(iu)
          write(iu, '(6x,a)') '</nn>'
        end if
        write(iu, '(6x,a)') '<mm>'
        do k = 1, elg%nel; do i = 1, FEDB(tp)%lnv; call feed(iu, string(elg%mm(i,k))); end do; end do; call empty(iu)
        write(iu, '(6x,a)') '</mm>'
        write(iu, '(6x,a)') '<ref>'
        do k = 1, elg%nel; call feed(iu, string(elg%ref(k))); end do; call empty(iu)
        write(iu, '(6x,a)') '</ref>'
      endif
      if(allocated(elg%fi)) then
        do ifi=1, size(elg%fi,1)
          if(.not. allocated(elg%fi(ifi)%val)) cycle
          write(iu, '(6x,a)') '<field name="'//trim(elg%fi(ifi)%name)//&
                             & '" ncomp="'//trim(string(size(elg%fi(ifi)%val,1)))//&
                             & '" nshot="'//trim(string(size(elg%fi(ifi)%param,1)))//'">'
          write(iu, '(6x,a)') '</field>'
        enddo
      endif
      write(iu, '(4x,a)') '</element_group>'
    end associate
  end do
  write(iu, '(2x,a)') '</piece>'
end do
write(iu, '(a)') '</pmh>'
if(wv) close(iu)
end subroutine

!-----------------------------------------------------------------------
! pmh2mfm: convert a PMH structure into a MFM one
!
! pmh is deallocated while MFM variables are being allocated
!-----------------------------------------------------------------------
subroutine pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
type(pmh_mesh), intent(inout) :: pmh
integer, intent(inout) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
integer, allocatable   :: nn(:,:), mm(:,:), nrv(:,:), nra(:,:), nrc(:,:), nsd(:)
real(real64), allocatable :: z(:,:)

integer :: i, ipp, ip, ig, pos, k, prev_nel, n, j, type_by_tdim(0:3), tmp_2d(2), tmp_3d(3), &
prev_max_tdim, res, max_tdim, valid_fe(13)
integer, allocatable :: piece2save(:), ref(:,:), tmp_vf(:), nel_piece(:), nnod_piece(:), nver_piece(:)
logical :: ft(2) = [.false.,.true.]
character(maxpath) :: str, cad

!valid elements types to save a MFM mesh (all but wedges)
valid_fe = [check_fe(.true.,   1, 1,  0, 0), & !Node
            check_fe(.true.,   2, 2,  1, 0), & !Edge, Lagrange P1
            check_fe(.false.,  3, 2,  1, 0), & !Edge, Lagrange P2
            check_fe(.true.,   3, 3,  3, 0), & !Triangle, Lagrange P1
            check_fe(.false.,  6, 3,  3, 0), & !Triangle, Lagrange P2
            check_fe(.false.,  3, 3,  3, 0), & !Triangle, Raviart-Thomas (edge)
            check_fe(.true.,   4, 4,  4, 0), & !Quadrangle, Lagrange P1
            check_fe(.true.,   4, 4,  6, 4), & !Tetrahedron, Lagrange P1
            check_fe(.false., 10, 4,  6, 4), & !Tetrahedron, Lagrange P2
            check_fe(.false.,  4, 4,  6, 4), & !Tetrahedron, Raviart-Thomas (face)
            check_fe(.false.,  6, 4,  6, 4), & !Tetrahedron, Nedelec (edge)
            check_fe(.false., 20, 4,  6, 4), & !Tetrahedron, Nedelec (edge) order 2
            check_fe(.true.,   8, 8, 12, 6)]   !Hexahedron, Lagrange P1

!check piece(s) to be saved
if (is_arg('-p')) then !save a single piece, indicated after -p
  str = get_post_arg('-p')
  call alloc(piece2save, 1)
  piece2save(1) = int(str)
else !save all pieces
  call alloc(piece2save, size(pmh%pc,1))
  piece2save = [(i, i=1, size(pmh%pc,1))]
end if
if (is_arg('-glue')) then
  call info('(module_pmh/pmh2mfm) option -glue not implemented yet')
end if

!testing and calculation of max_tdim
type_by_tdim  = 0 !store the type of element for each topological dimension
prev_max_tdim = 0 !store the maximal topological dimension
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_pmh/pmh2mfm) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (find_first(valid_fe, tp) == 0) then
        call info('(module_pmh/pmh2mfm) element type '//trim(FEDB(tp)%desc)//' found; those elements cannot be saved'//&
        &' in MFM format and they will be discarded')
        cycle
      end if
      !check whether there is only one type of element for each topological dimension
      if (type_by_tdim( FEDB(tp)%tdim ) == 0) then
        type_by_tdim( FEDB(tp)%tdim ) = tp
      elseif (type_by_tdim( FEDB(tp)%tdim ) /= tp) then
        call error('(module_pmh/pmh2mfm) more that one type of element is defined for the same topological dimension: '//&
        &string(type_by_tdim(FEDB(tp)%tdim))//', '//string(tp)//'; unable to convert to MFM')
      end if
    end associate
  end do
  max_tdim = 0
  do i = 1, 3
    if (type_by_tdim(i) > 0) max_tdim  = i
  end do
  if (prev_max_tdim == 0) then
    prev_max_tdim = max_tdim
  elseif (prev_max_tdim /= max_tdim) then
    call error('(module_pmh/pmh2mfm) there are pieces with different maximal topological dimension; unable to convert to MFM')
  end if
end do

!testing and calculation of dim
dim = -1
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  !check whether there is only one coordinates dimension for all pieces
  if (dim == -1) then
    dim = pmh%pc(ip)%dim
  elseif (pmh%pc(ip)%dim /= dim) then
    call error('(module_pmh/pmh2mfm) different coordinates dimensions in different pieces: '//&
    &string(pmh%pc(ip)%dim)//', '//string(dim)//'; unable to convert to MFM')
  end if
end do

!save MFM variables (lnn, lnv, lne, lnf) for selected pieces
!(previous testing ensure that lnn, lnv, etc., are the same for each piece)
lnn  = FEDB(type_by_tdim(max_tdim))%lnn
lnv  = FEDB(type_by_tdim(max_tdim))%lnv
lne  = FEDB(type_by_tdim(max_tdim))%lne
lnf  = FEDB(type_by_tdim(max_tdim))%lnf

!save MFM variables (nnod, nver, nel) for selected pieces
if (allocated(nel_piece)) deallocate(nel_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(nel_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate variable nel_piece: '//trim(cad))
if (allocated(nnod_piece)) deallocate(nnod_piece); !nnod_piece(ipp): global numbering for the last node  of piece #ipp
allocate(nnod_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate variable nnod_piece: '//trim(cad))
if (allocated(nver_piece)) deallocate(nver_piece); !nver_piece(ipp): global numbering for the last vertex of piece #ipp
allocate(nver_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate variable nver_piece: '//trim(cad))
nel_piece(0) = 0; nnod_piece(0) = 0; nver_piece(0) = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  nnod_piece(ipp) = nnod_piece(ipp-1) + pmh%pc(ip)%nnod
  nver_piece(ipp) = nver_piece(ipp-1) + pmh%pc(ip)%nver
  nel_piece(ipp)  =  nel_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (find_first(valid_fe, elg%type) == 0) cycle
      if (FEDB( elg%type )%tdim == max_tdim) then
         nel_piece(ipp) =  nel_piece(ipp) + elg%nel
      end if
    end associate
  end do
end do
nel  =  nel_piece(size(piece2save,1))
nver = nver_piece(size(piece2save,1))
nnod = nnod_piece(size(piece2save,1))

!save mm: concatenate vertex numbering of maximal topological dimension
call alloc(mm, lnv, nel)
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  prev_nel = nel_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (find_first(valid_fe, elg%type) == 0) cycle
      if (FEDB( elg%type )%tdim == max_tdim) then
        mm(1:lnv, prev_nel+1:prev_nel+elg%nel) = nver_piece(ipp-1) + elg%mm
        call dealloc(elg%mm)
        prev_nel = prev_nel + elg%nel
      end if
    end associate
  end do
end do

!nsd: concatenate references of maximal topological dimension
call alloc(nsd, nel); nsd = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  prev_nel = nel_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (find_first(valid_fe, elg%type) == 0) cycle
      if (FEDB( elg%type )%tdim == max_tdim) then
        nsd(prev_nel+1:prev_nel+elg%nel) = elg%ref
        call dealloc(elg%ref)
        prev_nel = prev_nel + elg%nel
      end if
    end associate
  end do
end do

!nn: concatenate node numbering of maximal topological dimension
if (.not. FEDB(type_by_tdim(max_tdim))%nver_eq_nnod) then
  call alloc(nn, lnn, nel)
  do ipp = 1, size(piece2save,1)
    ip = piece2save(ipp)
    prev_nel = nel_piece(ipp-1)
    do ig = 1, size(pmh%pc(ip)%el, 1)
      associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
        if (find_first(valid_fe, elg%type) == 0) cycle
        if (FEDB( elg%type )%tdim == max_tdim) then
          nn(1:lnn, prev_nel+1:prev_nel+elg%nel) = nnod_piece(ipp-1) + elg%nn
          call dealloc(elg%nn)
          prev_nel = prev_nel + elg%nel
        end if
      end associate
    end do
  end do
end if

!nrv: visit element groups of tdim = 0 to set vertex references
if (max_tdim > 0) then
  call alloc(nrv, lnv, nel); nrv = 0
  do ipp = 1, size(piece2save,1)
    ip = piece2save(ipp)
    !STEP 1: create ref to collect vertices and references stored in PMH, orderly
    n = 0
    do ig = 1, size(pmh%pc(ip)%el, 1)
      associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
        if (find_first(valid_fe, elg%type) == 0) cycle
        if (FEDB( elg%type )%tdim == 0) then
          do k = 1, elg%nel
            tmp_2d = [nver_piece(ipp-1) + elg%mm(1,k), elg%ref(k)]
            call insert_sorted(1, ref, tmp_2d, used=n, fit=ft)
          end do
          call dealloc(elg%ref)
        end if
      end associate
    end do
    if (allocated(ref)) then
      call reduce(ref, n, 2)
      !STEP2: for every vertex in mm, check whether it is in ref
      do k = nel_piece(ipp-1)+1, nel_piece(ipp)
        do j = 1, FEDB(type_by_tdim(max_tdim))%lnv
          pos = find_first(ref(:,1), mm(j,k))
          if (pos > 0) nrv(j,k) = ref(pos, 2)
        end do
      end do
      call dealloc(ref)
    end if
  end do
end if

!nra: visit element groups of tdim = 1 to set edge references
if (max_tdim > 1) then
  call alloc(nra, lne, nel); nra = 0
  do ipp = 1, size(piece2save,1)
    ip = piece2save(ipp)
    !STEP 1: create ref to collect edge vertices and references stored in PMH, orderly
    n = 0
    do ig = 1, size(pmh%pc(ip)%el, 1)
      associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
        if (find_first(valid_fe, elg%type) == 0) cycle
        if (FEDB( elg%type )%tdim == 1) then
          do k = 1, elg%nel
            tmp_3d = [sort(nver_piece(ipp-1) + elg%mm(1:2,k)), elg%ref(k)]
            call insert_sorted(1, ref, tmp_3d, used=n, fit=ft)
          end do
          call dealloc(elg%ref)
        end if
      end associate
    end do
    if (allocated(ref)) then
      call reduce(ref, n, 3)
      !STEP2: for every edge in mm, check whether it is in ref
      do k = nel_piece(ipp-1)+1, nel_piece(ipp)
        do j = 1, FEDB(type_by_tdim(max_tdim))%lne
          tmp_2d = sort(mm(FEDB(type_by_tdim(max_tdim))%edge(:,j),k))
          pos = find_sorted(1, ref, tmp_2d, n)
          if (pos > 0) nra(j,k) = ref(pos, 3)
        end do
      end do
      call dealloc(ref)
    end if
  end do
end if
!nrc: visit element groups of tdim = 2 to set edge references
if (max_tdim > 2) then
  call alloc(nrc, lnf, nel); nrc = 0
  do ipp = 1, size(piece2save,1)
    ip = piece2save(ipp)
    !STEP 1: create ref2 to collect face vertices and references stored in PMH, orderly
    n = 0
    associate(v_f => FEDB(FEDB(type_by_tdim(max_tdim))%f_type)%lnv) !v_f: vertices per face in the max. tdim element
      call alloc(tmp_vf, v_f+1)
      do ig = 1, size(pmh%pc(ip)%el, 1)
        associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
          if (find_first(valid_fe, elg%type) == 0) cycle
          if (FEDB( elg%type )%tdim == 2 .and. FEDB(elg%type)%lnv == v_f) then
            do k = 1, elg%nel
              tmp_vf = [sort(nver_piece(ipp-1) + elg%mm(1:v_f, k)), elg%ref(k)]
              call insert_sorted(1, ref, tmp_vf, used=n, fit=ft)
            end do
            call dealloc(elg%ref)
          end if
        end associate
      end do
      if (allocated(ref)) then
        call reduce(ref, n, v_f+1)
        !STEP2: for every face in mm, check whether it is in ref
        call alloc(tmp_vf, v_f)
        do k = nel_piece(ipp-1)+1, nel_piece(ipp)
          do j = 1, FEDB(type_by_tdim(max_tdim))%lnf
            tmp_vf = sort(mm(FEDB(type_by_tdim(max_tdim))%face(1:v_f, j), k))
            pos = find_sorted(1, ref, tmp_vf, n)
            if (pos > 0) nrc(j,k) = ref(pos, v_f+1)
          end do
        end do
        call dealloc(ref)
      end if
    end associate
  end do
end if

!z: save vertex coordinates
call alloc(z, dim, nver)
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  z(1:dim, nver_piece(ipp-1)+1:nver_piece(ipp)) = pmh%pc(ip)%z
  deallocate(pmh%pc(ip)%z)
end do
end subroutine

!-----------------------------------------------------------------------
! mfm2pmh: convert a MFM mesh into a PMH structure
!
! mfm is deallocated while PMH variables are being allocated
!-----------------------------------------------------------------------
subroutine mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh, deallocation)
type(pmh_mesh), intent(inout) :: pmh
integer, intent(inout) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
integer, allocatable   :: nn(:,:), mm(:,:), nrv(:,:), nra(:,:), nrc(:,:), nsd(:)
real(real64), allocatable :: z(:,:)
logical, optional, intent(in) :: deallocation
integer :: res, tp, nelg, iel, ic, j, k
character(maxpath) :: cad
logical :: deall

!whether or not to deallocate the MFM variables (.true. by default)
deall = .true.
if (present(deallocation)) deall = deallocation

!allocate piece, pmh%pc
if (allocated(pmh%pc)) then
  if (size(pmh%pc,1) /= 1) then
    deallocate(pmh%pc, stat = res, errmsg = cad)
    if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to deallocate piece: '//trim(cad))
    allocate(pmh%pc(1), stat = res, errmsg = cad)
    if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate piece: '//trim(cad))
  end if
else
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate piece: '//trim(cad))
end if

!save PMH variables (nnod, nver, dim, z)
pmh%pc(1)%nnod = nnod
pmh%pc(1)%nver = nver
pmh%pc(1)%dim  = dim
if (deall) then
  call move_alloc(from=z, to=pmh%pc(1)%z)
else
  call alloc(pmh%pc(1)%z, size(z,1), size(z,2))
  pmh%pc(1)%z = z
end if

!nelg: calculate the number of element groups to create
tp = check_fe(nver==nnod, lnn, lnv, lne, lnf)
if (tp == 0) then
  call error('(module_pmh/pmh2mfm) Unable to process a MFM mesh with nver: '//trim(string(nver))//', nnod: '//&
  &trim(string(nnod))//', lnn: '//trim(string(lnn))//', lnv: '//trim(string(lnv))//', lne:'//trim(string(lne))//', lnf: '//&
  &trim(string(lnf))//'. Please check the valid finite elements in the structure FEDB of source/pmh/module_fe_database_pmh.f90')
end if
nelg = FEDB(tp)%tdim + 1

if (FEDB(tp)%tdim > 2) then
  if (maxval(nrc)==0) nelg = nelg-1
end if
if (FEDB(tp)%tdim > 1) then
  if (maxval(nra)==0) nelg = nelg-1
end if
if (maxval(nrv)==0) nelg = nelg-1

!allocate element groups, pmh%pc(1)%el
if (allocated(pmh%pc(1)%el)) then
  if (size(pmh%pc(1)%el,1) /= nelg) then
    deallocate(pmh%pc(1)%el, stat = res, errmsg = cad)
    if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to deallocate element groups: '//trim(cad))
    allocate(pmh%pc(1)%el(nelg), stat = res, errmsg = cad)
    if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate element groups: '//trim(cad))
  end if
else
  allocate(pmh%pc(1)%el(nelg), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_pmh/pmh2mfm) Unable to allocate element groups: '//trim(cad))
end if

!highest tdim: mm, nn and nsd
associate (elg => pmh%pc(1)%el(nelg)) !elg: current element group
  elg%nel = nel
  elg%type = tp
  if (deall) then
    call move_alloc(from=nn,  to=pmh%pc(1)%el(nelg)%nn)
    call move_alloc(from=mm,  to=pmh%pc(1)%el(nelg)%mm)
    call move_alloc(from=nsd, to=pmh%pc(1)%el(nelg)%ref)
  else
    call alloc(pmh%pc(1)%el(nelg)%nn, size(nn,1), size(nn,2)); pmh%pc(1)%el(nelg)%nn  = nn
    call alloc(pmh%pc(1)%el(nelg)%mm, size(mm,1), size(mm,2)); pmh%pc(1)%el(nelg)%mm  = mm
    call alloc(pmh%pc(1)%el(nelg)%ref, size(nsd,1));           pmh%pc(1)%el(nelg)%ref = nsd
  end if
end associate
iel = nelg - 1

!group for faces
if (FEDB(tp)%tdim > 2) then
  if (maxval(nrc) > 0) then
    associate (melg => pmh%pc(1)%el(nelg), elg => pmh%pc(1)%el(iel)) !melg: max. tdim group, elg: current group
      elg%nel = count(nrc > 0)
      elg%type = FEDB(melg%type)%f_type
      !in faces, nn is only relevant for Lagrange P2
      if (FEDB(elg%type)%lnn == FEDB(elg%type)%lnv + FEDB(elg%type)%lne) &
      call alloc(elg%nn, FEDB(elg%type)%lnn, elg%nel)
      call alloc(elg%mm, FEDB(elg%type)%lnv, elg%nel)
      call alloc(elg%ref, elg%nel)
      ic = 1
      do k = 1, melg%nel
        do j = 1, FEDB(melg%type)%lnf
          if (nrc(j,k) /= 0) then
            if (FEDB(elg%type)%lnn == FEDB(elg%type)%lnv + FEDB(elg%type)%lne) &
            elg%nn(:, ic) = melg%nn(FEDB(melg%type)%nface(1:FEDB(elg%type)%lnn, j), k)
            elg%mm(:, ic) = melg%mm(FEDB(melg%type)%face( 1:FEDB(elg%type)%lnv, j), k)
            elg%ref(ic)   = nrc(j,k)
            ic = ic + 1
          end if
        end do
      end do
    end associate
    iel = iel - 1
  end if
  if (deall) call dealloc(nrc)
end if

!group for edges
if (FEDB(tp)%tdim > 1) then
  if (maxval(nra) > 0) then
    associate (melg => pmh%pc(1)%el(nelg), elg => pmh%pc(1)%el(iel)) !melg: max. tdim group, elg: current group
      elg%nel = count(nra > 0)
      elg%type = FEDB(melg%type)%e_type
      !in edges, nn is only relevant for Lagrange P2
      if (FEDB(elg%type)%lnn == FEDB(elg%type)%lnv + FEDB(elg%type)%lne) &
      call alloc(elg%nn, FEDB(elg%type)%lnn, elg%nel)
      call alloc(elg%mm, FEDB(elg%type)%lnv, elg%nel)
      call alloc(elg%ref, elg%nel)
      ic = 1
      do k = 1, melg%nel
        do j = 1, FEDB(melg%type)%lne
          if (nra(j,k) /= 0) then
            if (FEDB(elg%type)%lnn == FEDB(elg%type)%lnv + FEDB(elg%type)%lne) &
            elg%nn(:, ic) = melg%nn(FEDB(melg%type)%nedge(1:FEDB(elg%type)%lnn, j), k)
            elg%mm(:, ic) = melg%mm(FEDB(melg%type)%edge( 1:FEDB(elg%type)%lnv, j), k)
            elg%ref(ic)   = nra(j,k)
            ic = ic + 1
          end if
        end do
      end do
    end associate
    iel = iel - 1
  end if
  if (deall) call dealloc(nra)
end if

!group for vertices
if (FEDB(tp)%tdim > 0) then
  if (maxval(nrv) > 0) then
    associate (melg => pmh%pc(1)%el(nelg), elg => pmh%pc(1)%el(iel)) !melg: max. tdim group, elg: current group
      elg%nel = count(nrv > 0)
      elg%type = FEDB(melg%type)%v_type
      !in vertices, nn is not relevant
      call alloc(elg%mm, FEDB(elg%type)%lnv, elg%nel)
      call alloc(elg%ref, elg%nel)
      ic = 1
      do k = 1, melg%nel
        do j = 1, FEDB(melg%type)%lnv
          if (nrv(j,k) /= 0) then
            elg%mm(1, ic) = melg%mm(j,k)
            elg%ref(ic)   = nrv(j,k)
            ic = ic + 1
          end if
        end do
      end do
    end associate
  end if
  if (deall) call dealloc(nrv)
end if

end subroutine

!-----------------------------------------------------------------------
! build_vertices: build vertex connectivities and coordinates from node information
!
! vertex built is only done for Lagrange P1 and P2 elements
! it is assumed that vertices are the first nodes, ensured by reorder_nodes()
! for P2, it is assumed that global numbering as vertex is less or equal than global numbering as node
!-----------------------------------------------------------------------
subroutine build_vertices(pmh)
type(pmh_mesh), intent(inout) :: pmh
integer, allocatable :: vert2node(:), node2vert(:)
integer :: nv2d, pos, maxv, i, j, k, ig, ip
logical :: all_are_P1
logical, allocatable :: unalloc_mm_P2(:)

!reorder_nodes: reorder nodes and/or vertices of Lagrange P2 elements to have vertices before mid-points
call reorder_nodes_P2(pmh)

!create mm (and reconstruct z) if necessary
do ip = 1, size(pmh%pc,1)
  if (.not. allocated(pmh%pc(ip)%z)) call error('(module_pmh/build_vertices) z is not allocated: piece '//trim(string(ip))//&
  &'; unable to build vertices')
  if (pmh%pc(ip)%nver /= size(pmh%pc(ip)%z,2) .and. pmh%pc(ip)%nnod /= size(pmh%pc(ip)%z,2)) &
  call error('(module_pmh/build_vertices) z has an incorrect number of columns: '//trim(string(size(pmh%pc(ip)%z,2)))//&
  &' while nver is '//trim(string(pmh%pc(ip)%nver))//' and nnod is '//trim(string(pmh%pc(ip)%nnod))//': piece '//trim(string(ip)))
  all_are_P1 = .true.
  call alloc(unalloc_mm_P2, size(pmh%pc(ip)%el,1))
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
      if (FEDB(tp)%nver_eq_nnod) then
        !Lagrange P1: create mm and erase nn
        if (.not. allocated(elg%mm)) then
          if (.not.allocated(elg%nn)) call error('(module_pmh/build_vertices) neither mm nor nn are allocated: piece '//&
          &trim(string(ip))//', group '//trim(string(ig))//'; unable to build vertices')
          call move_alloc(from=elg%nn, to=elg%mm)
        end if
        call dealloc(elg%nn)
      elseif (FEDB(tp)%lnn == FEDB(tp)%lnv + FEDB(tp)%lne) then
        !Lagrange P2: if mm is not defined, it will created later with vert2node
        if (.not. allocated(elg%mm)) then
          if (.not.allocated(elg%nn)) call error('(module_pmh/build_vertices) neither mm nor nn are allocated: piece '//&
          &trim(string(ip))//', group '//trim(string(ig))//'; unable to build vertices')
          unalloc_mm_P2(ig) = .true.
        end if
        all_are_P1 = .false.
      else
        !Other types: if mm is not defined, it cannot be created from nn
        if (.not. allocated(elg%mm)) call error('(module_pmh/build_vertices) mm is not allocated and there are elements of type '//&
        &trim(FEDB(tp)%desc)//': piece '//trim(string(ip))//', group '//trim(string(ig))//'; unable to build vertices')
        all_are_P1 = .false.
      end if
    end associate
  end do
  if (any(unalloc_mm_P2) .or. (.not. all_are_P1 .and. pmh%pc(ip)%nnod == size(pmh%pc(ip)%z,2))) then
    !there are some unallocated mm for Lagrange P2 elements or z must be reconstructed: construct global vertex numbering
    nv2d = 0
    do ig = 1, size(pmh%pc(ip)%el,1)
      associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnv
            if (allocated(elg%mm)) then
              !elements have mm already constructed: save them in vert2node
              pos = bsearch(vert2node, elg%mm(i,k), nv2d)
              if (pos < 0) then
                call insert(vert2node, elg%mm(i,k), -pos, nv2d, fit=.false.)
                nv2d = nv2d + 1
              end if
            else
              !only nn is allocated; previous checkings ensure that it must be a Lagrange P2 element
              if (elg%nn(i,k) == 0) call error('(module_pmh/build_vertices) node numbering is zero: piece '//trim(string(ip))//&
              &', group '//trim(string(ig))//', node '//trim(string(i))//', element '//trim(string(k))//'; unable to build vertex')
              pos = bsearch(vert2node, elg%nn(i,k), nv2d)
              if (pos < 0) then
                call insert(vert2node, elg%nn(i,k), -pos, nv2d, fit=.false.)
                nv2d = nv2d + 1
              end if
            end if
          end do
        end do
      end associate
    end do
    call reduce(vert2node, nv2d)
    !nver: total number of vertices
    pmh%pc(ip)%nver = nv2d
    !node2vert: given the global numbering of vertices (numbered as nodes), returns global numbering as vertices
    maxv = 0
    do i = 1, size(vert2node, 1)
      maxv = max(maxv, vert2node(i))
    enddo
    call alloc(node2vert, maxv)
    do i = 1, size(vert2node, 1)
      if (vert2node(i) == 0) cycle
      node2vert(vert2node(i)) = i
    enddo
    if (any(unalloc_mm_P2)) then
      !only if there are some unallocated mm for Lagrange P2 elements, we use global vertex numbering to create mm
      do ig = 1, size(pmh%pc(ip)%el,1)
        if (.not. unalloc_mm_P2(ig)) cycle !only for groups of Lagrange P2 elements where mm is not allocated
          associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
            call alloc(elg%mm, FEDB(tp)%lnv, elg%nel)
            do k = 1, elg%nel
              do i = 1, FEDB(tp)%lnv
              elg%mm(i,k) = node2vert(elg%nn(i,k))
            end do
          end do
        end associate
      end do
    end if
    !vertex coordinates
    associate(pi => pmh%pc(ip))
      if (pi%nver == size(pi%z,2)) then
        call info('(module_pmh/build_vertices) z has the correct number of columns although some mm for Lagrange P2 '//&
        &'elements where initially unallocated')
      elseif (pi%nnod == size(pi%z,2)) then !z stores node coord.: assume that i <= vert2node(i); thus z can be overwritten
        do j = 1, pi%nver
          pi%z(1:pi%dim, j) = pi%z(1:pi%dim, vert2node(j))
        end do
        call reduce(pi%z, pi%dim, pi%nver)
      else
        call error('(module_pmh/build_vertices (2)) z has an incorrect number of columns: '//trim(string(size(pmh%pc(ip)%z,2)))//&
        &' while nver is '//trim(string(pmh%pc(ip)%nver))//' and nnod is '//trim(string(pmh%pc(ip)%nnod))//': piece '//&
        &trim(string(ip)))
      end if
    end associate
    call dealloc(vert2node)
    call dealloc(node2vert)
  end if
  !at the end, z must store only vertex coordinates
  if (pmh%pc(ip)%nver == 0) then
    if (all_are_P1) then
      pmh%pc(ip)%nver = pmh%pc(ip)%nnod
    else
      call error('(module_pmh/build_vertices) nver is still zero: piece '//trim(string(ip)))
    end if
  end if
  if (pmh%pc(ip)%nver /= size(pmh%pc(ip)%z,2)) &
  call error('(module_pmh/build_vertices (3)) z has an incorrect number of columns: '//trim(string(size(pmh%pc(ip)%z,2)))//&
  &' while nver is '//trim(string(pmh%pc(ip)%nver))//' and nnod is '//trim(string(pmh%pc(ip)%nnod))//': piece '//trim(string(ip)))
end do

!positive_jacobian: ensure positive jacobian
call positive_jacobian(pmh)
end subroutine

!-----------------------------------------------------------------------
! build_node_coordinates: build node coordinates from vertex information
!
! this procedure must be called for each piece of a well-constructed PMH
! When all grous are Lagrange P1, variable all_P1 is .true. and znod is not created
! thus, the way of calling this procedure is:
!   do ip = 1, size(pmh%pc, 1)
!     call build_node_coordinates(pmh%pc(ip), ip, all_P1, znod)
!     if (.not. all_P1) then
!        ... work with znod ...
!     end if
!   end do
!-----------------------------------------------------------------------
subroutine build_node_coordinates(pc, ip, all_P1, znod)
type(piece),               intent(in)  :: pc
integer,                   intent(in)  :: ip
logical,                   intent(out) :: all_P1
real(real64), allocatable, intent(out) :: znod(:,:)
integer :: ig, i, j, k

!determine whether all elements are Lagrange P1
all_P1 = .true.
do ig = 1, size(pc%el,1)
  if (.not. FEDB(pc%el(ig)%type)%nver_eq_nnod) then
    all_P1 = .false.
    exit
  end if
end do
if (.not. all_P1) then
  !construct znod
  call alloc(znod, pc%dim, pc%nnod)
  do ig = 1, size(pc%el,1)
    associate(elg => pc%el(ig), tp => pc%el(ig)%type)
      if (FEDB(tp)%nver_eq_nnod) then
        !************************************* Lagrange P1 **************************************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnv
            znod(:, elg%mm(i,k)) = pc%z(:, elg%mm(i,k))
          end do
        end do
      elseif (FEDB(tp)%lnn == FEDB(tp)%lnv + FEDB(tp)%lne) then
        !************************************* Lagrange P2 **************************************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnv
            znod(:, elg%nn(i,k)) = pc%z(:, elg%mm(i,k))
          end do
          do i = 1, FEDB(tp)%lne
            znod(:, elg%nn(i+FEDB(tp)%lnv,k)) = (pc%z(:,elg%mm(FEDB(tp)%edge(1,i),k)) + pc%z(:,elg%mm(FEDB(tp)%edge(2,i),k)))/2
          end do
        end do
      elseif (tp == check_fe(.false.,  3, 3,  3, 0) .or. &
              tp == check_fe(.false.,  6, 4,  6, 4)) then
        !***** Triangle, Raviart-Thomas (edge) OR Tetrahedron, Nedelec (edge) ******************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lne
            znod(:, elg%nn(i,k)) = ( pc%z(:, elg%mm(FEDB(tp)%edge(1,i),k)) + pc%z(:, elg%mm(FEDB(tp)%edge(2,i),k)) )/2
          end do
        end do
      elseif (tp == check_fe(.false.,  4, 4,  6, 4)) then
        !************************************* Tetrahedron, Raviart-Thomas (face) ***************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnf
            do j = 1, pc%dim
              znod(j, elg%nn(i,k)) = sum(pc%z(j, elg%mm(FEDB(tp)%face(:,i),k))) / FEDB(FEDB(tp)%f_type)%lnv
            end do
          end do
        end do
      else
        call info('(module_pmh/build_node_coordinates) build node coordinates for element type '//trim(string(FEDB(tp)%desc))//&
        &' is not implemented: piece '//trim(string(ip))//', group '//trim(string(ig))//'; variable znod is not created')
      end if
    end associate
  end do
end if
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! reorder_nodes_P2: ensure that mid-points in Lagrange P2 appear at the
! end of conectivity arrays
!
! it must be called at the beginning of build_vertices
!-----------------------------------------------------------------------
subroutine reorder_nodes_P2(pmh)
type(pmh_mesh), intent(inout) :: pmh
integer :: ip, ig, k
integer, allocatable :: inew(:)
character(maxpath) :: reorder_type

!check what type of reorder must be applied
if (is_arg('-r')) then
  select case(get_post_arg('-r'))
  case('hard', 'soft', 'sandwich') !recognized options
    reorder_type = get_post_arg('-r')
  case default
    call error('(module_pmh/reorder_nodes) option -r not recognized: '//trim(get_post_arg('-r'))//&
    &'; use ''feconv -h'' to see available options')
  end select
else
  reorder_type = 'hard'
end if
do ip = 1, size(pmh%pc,1)
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type, z => pmh%pc(ip)%z)
      if (FEDB(tp)%nver_eq_nnod) then !Lagrange P1 elements, do nothing
        continue
      elseif (tp == check_fe(.false.,  3, 2,  1, 0)) then
        !************************************* Edge, Lagrange P2 ********************************
        !reorder nn such that mid-points appear at the end
        if (.not.allocated(elg%nn)) call error('(module_pmh/reorder_nodes) nn not allocated for an Lagrange P2 edge mesh: piece '//&
        &trim(string(ip))//', group '//trim(string(ig)))
        select case(reorder_type)
        case ('hard') !check every element
          do k = 1, elg%nel
            if     (maxval(abs((z(:,elg%nn(1,k))+z(:,elg%nn(2,k)))/2-z(:,elg%nn(3,k)))) < 1e3*pmh%ztol) then
              elg%nn([1, 2, 3],k) = elg%nn(:,k)
            elseif (maxval(abs((z(:,elg%nn(1,k))+z(:,elg%nn(3,k)))/2-z(:,elg%nn(2,k)))) < 1e3*pmh%ztol) then
              elg%nn([1, 3, 2],k) = elg%nn(:,k)
            elseif (maxval(abs((z(:,elg%nn(2,k))+z(:,elg%nn(3,k)))/2-z(:,elg%nn(1,k)))) < 1e3*pmh%ztol) then
              elg%nn([2, 3, 1],k) = elg%nn(:,k)
            else
              call info('(module_pmh/reorder_nodes) the nodes of an Lagrange P2 edge do not lie on a straight line: piece '//&
              &trim(string(ip))//', group '//trim(string(ig))//', element '//trim(string(k))//'; reordering is not guaranteed.')
            end if
          end do
        case('soft') !check the first element, asume the rest have the same behavior
          call alloc(inew, FEDB(tp)%lnn)
          k = 1
          if     (maxval(abs((z(:,elg%nn(1,k))+z(:,elg%nn(2,k)))/2-z(:,elg%nn(3,k)))) < 1e3*pmh%ztol) then
            inew = [1, 2, 3]
          elseif (maxval(abs((z(:,elg%nn(1,k))+z(:,elg%nn(3,k)))/2-z(:,elg%nn(2,k)))) < 1e3*pmh%ztol) then
            inew = [1, 3, 2]
          elseif (maxval(abs((z(:,elg%nn(2,k))+z(:,elg%nn(3,k)))/2-z(:,elg%nn(1,k)))) < 1e3*pmh%ztol) then
            inew = [2, 3, 1]
          else
            call info('(module_pmh/reorder_nodes) the nodes of an Lagrange P2 edge do not lie on a straight line: piece '//&
            &trim(string(ip))//', group '//trim(string(ig))//', element '//trim(string(k))//'; reordering is not guaranteed.')
          end if
          do k = 1, elg%nel
            elg%nn(inew(:),k) = elg%nn(:,k)
          end do
        case('sandwich') !node order is the one prescribed by SALOME: vert. and mid-points sandwiched
          do k = 1, elg%nel
            elg%nn([1,3,2],k) = elg%nn(:,k)
          end do
        end select
      elseif (tp == check_fe(.false.,  6, 3, 3, 0) .or. tp == check_fe(.false.,  8, 4, 4, 0)) then
        !************************************* Triangle and quadrangle, Lagrange P2 ****************************
        if (.not.allocated(elg%nn)) &
          & call error('(module_pmh/reorder_nodes) nn not allocated for an Lagrange P2 surface mesh: piece '//&
        &trim(string(ip))//', group '//trim(string(ig)))
        select case(reorder_type)
        case('hard') !check every element
          !check whether mid-points appear at the end of nn
          call alloc(inew, FEDB(tp)%lnn)
          do k = 1, elg%nel
            call reorder_nodes_element_P2(ip, ig, elg, tp, z, pmh%ztol, k, inew)
            elg%nn(:,k) = elg%nn(inew(:),k)
          end do
        case('soft') !check the first element, asume the rest have the same behavior
          !check whether mid-points appear at the end of nn
          call alloc(inew, FEDB(tp)%lnn)
          call reorder_nodes_element_P2(ip, ig, elg, tp, z, pmh%ztol, 1, inew)
          do k = 1, elg%nel
            elg%nn(:,k) = elg%nn(inew(:),k)
          end do
        case('sandwich') !node order is the one prescribed by SALOME for isoparam. P2 triangles: vert. and mid-points sandwiched
          do k = 1, elg%nel
            if(tp == check_fe(.false.,  6, 3, 3, 0)) then
              elg%nn([1,4,2,5,3,6],k) = elg%nn(:,k)                !Tria P2
            elseif(tp == check_fe(.false.,  8, 4, 4, 0)) then
              elg%nn([1,5,2,6,3,7,4,8],k) = elg%nn(:,k)            !Quad P2
            endif
          end do
        end select
      elseif (tp == check_fe(.false., 10, 4,  6, 4) .or. tp == check_fe(.false., 20, 8,  12, 6)) then
        !************************************* Tetrahedron and hexahedron, Lagrange P2 *************************
        if (.not.allocated(elg%nn)) &
          & call error('(module_pmh/reorder_nodes) nn not allocated for an Lagrange P2 volumic mesh: piece '//&
        &trim(string(ip))//', group '//trim(string(ig)))
        select case(reorder_type)
        case('hard') !check every element
          !check whether mid-points appear at the end of nn
          call alloc(inew, FEDB(tp)%lnn)
          do k = 1, elg%nel
            call reorder_nodes_element_P2(ip, ig, elg, tp, z, pmh%ztol, k, inew)
            elg%nn(:,k) = elg%nn(inew(:),k)
          end do
        case('soft') !check the first element, asume the rest have the same behavior
          !check whether mid-points appear at the end of nn
          call alloc(inew, FEDB(tp)%lnn)
          call reorder_nodes_element_P2(ip, ig, elg, tp, z, pmh%ztol, 1, inew)
          do k = 1, elg%nel
            elg%nn(:,k) = elg%nn(inew(:),k)
          end do
        case('sandwich') !node order is the one prescribed by SALOME for isoparam. P2 triangles: vert. and mid-points sandwiched
          do k = 1, elg%nel
            if(tp == check_fe(.false., 10, 4,  6, 4)) then
              elg%nn([1,5,2,6,3,7,8,9,10,4],k) = elg%nn(:,k)                                ! Tetra  P2
            elseif(tp == check_fe(.false., 20, 8,  12, 6)) then
              elg%nn([1,9,2,10,3,11,4,12,13,14,15,16,5,17,6,18,7,19,8,20],k) = elg%nn(:,k)  ! Hexa P2
            endif
          enddo
        end select
      else
        call info('(module_pmh/reorder_nodes) reordering of element type '//trim(string(FEDB(pmh%pc(ip)%el(ig)%type)%desc))//&
        &' is not implemented: piece '//trim(string(ip))//', group '//trim(string(ig))//'; node order remains unchanged')
      end if
    end associate
  end do
end do
end subroutine

!-----------------------------------------------------------------------
! reorder_nodes_element_P2: calculate inew to reorder nn (vertices first)
!-----------------------------------------------------------------------
subroutine reorder_nodes_element_P2(ip, ig, elg, tp, z, ztol, k, inew)
type(elgroup), intent(in)    :: elg
integer,       intent(in)    :: ip, ig, tp, k
real(real64),  intent(in)    :: z(:,:)
real(real64),  intent(in)    :: ztol
integer,       intent(inout) :: inew(:)
integer :: i, j, l, nv, newnn(size(inew))

newnn = 0
!search first the vertices, checking that they are not mid-points
nv = 0
NODES: do i = 1, FEDB(tp)%lnn !nodes
  do j = 1, FEDB(tp)%lnn !first possible vertex
    if (j /= i) then
      do l = j+1, FEDB(tp)%lnn !second possible vertex
        if (l /= i) then
          if ( maxval(abs((z(:,elg%nn(j,k))+z(:,elg%nn(l,k)))/2-z(:,elg%nn(i,k)))) < 1e3*ztol ) cycle NODES
        end if
      end do
    end if
  end do
  !elg%nn(i,k) is not a mid-point, so it is a vertex
  nv = nv + 1
  newnn(nv) = elg%nn(i,k)
  inew (nv) = i
  if (nv > FEDB(tp)%lnv)  call error('(module_pmh/reorder_nodes_element_P2) too many vertices were found in a Lagrange P2 '//&
  &'element: piece '//trim(string(ip))//', group '//trim(string(ig))//', element '//trim(string(k))//&
  &'; some edges can be singular or it can be an isoparametric element. Use ''feconv -h'' to see available options')
end do NODES
if (nv < FEDB(tp)%lnv)  call error('(module_pmh/reorder_nodes_element_P2) too few vertices were found in a Lagrange P2 '//&
&'element: piece '//trim(string(ip))//', group '//trim(string(ig))//', element '//trim(string(k))//&
&'; tolerance could be higher than the precision of the mesh. Use ''feconv -h'' to see available options')
!identify mid-points and save it in PMH order
do i = 1, FEDB(tp)%lnn !nodes
  if (find_first(newnn, elg%nn(i,k)) > 0) cycle !it is a vertex
  do j = 1, FEDB(tp)%lne
    if (maxval(abs((z(:,newnn(FEDB(tp)%edge(1,j))) + z(:,newnn(FEDB(tp)%edge(2,j))))/2 - z(:,elg%nn(i,k)))) < 1e3*ztol) then
      inew (FEDB(tp)%lnv + j) = i
      exit
    end if
  end do
end do
end subroutine

!-----------------------------------------------------------------------
! positive_jacobian: ensure positive jacobian in groups with maximal topological dimension
!
! it must be called at the end of build_vertices to guarantee that mm exists
!-----------------------------------------------------------------------
subroutine positive_jacobian(pmh)
type(pmh_mesh), intent(inout) :: pmh
integer :: ip, ig, i, k, max_tdim
integer, allocatable :: pos(:)
real(real64) ::  s, t, st(2)
logical :: QJac(4) ,check
real(real64) :: zz(2)=[0._real64, 0._real64],zo(2)=[0._real64, 1._real64],oz(2)=[1._real64, 0._real64]
!integer, allocatable :: tmp(:)

!check whether the application of possitive jacobian (by default, positive jacobian is run)
if (is_arg('-j')) then
  select case(get_post_arg('-j'))
  case('yes')
    continue
  case('no')
    call info('(module_pmh/positive_jacobian) detected option -j no: jacobian is not checked')
    return
  case default
    call error('(module_pmh/positive_jacobian) option -j not recognized: '//trim(get_post_arg('-j'))//&
    &'; use ''feconv -h'' to see available options')
  end select
end if

max_tdim = 0 !maximal topological dimension
do ip = 1, size(pmh%pc,1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
      if (max_tdim < FEDB(pmh%pc(ip)%el(ig)%type)%tdim) max_tdim = FEDB(pmh%pc(ip)%el(ig)%type)%tdim
  end do
end do

do ip = 1, size(pmh%pc,1)
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type, z => pmh%pc(ip)%z)
      if (FEDB(tp)%tdim < max_tdim) cycle !only for groups with maximal topological dimension
      if (tp == check_fe(.true.,   1, 1,  0, 0) .or. tp == check_fe(.true.,   2, 2,  1, 0) .or. &
          tp == check_fe(.false.,  3, 2,  1, 0)) cycle !Node, Edge Lagrange P1 or P2, do nothing
      if (tp == check_fe(.true.,   3, 3,  3, 0)) then
        !************************************* Triangle, Lagrange P1 ****************************
        do k = 1, elg%nel
          if ( (z(1,elg%mm(2,k)) - z(1,elg%mm(1,k))) * (z(2,elg%mm(3,k)) - z(2,elg%mm(1,k))) & !only x and y coordinates are used
             - (z(2,elg%mm(2,k)) - z(2,elg%mm(1,k))) * (z(1,elg%mm(3,k)) - z(1,elg%mm(1,k))) < 0 ) &
          call swap(elg%mm(2,k), elg%mm(3,k))
        end do
      elseif (tp == check_fe(.false.,  6, 3,  3, 0)) then
        !************************************* Triangle, Lagrange P2 ****************************
        do k = 1, elg%nel
          if ( (z(1,elg%mm(2,k)) - z(1,elg%mm(1,k))) * (z(2,elg%mm(3,k)) - z(2,elg%mm(1,k))) & !only x and y coordinates are used
             - (z(2,elg%mm(2,k)) - z(2,elg%mm(1,k))) * (z(1,elg%mm(3,k)) - z(1,elg%mm(1,k))) < 0 ) then
            call swap(elg%mm(2,k), elg%mm(3,k))
            call swap(elg%nn(2,k), elg%nn(3,k))
            call swap(elg%nn(4,k), elg%nn(6,k))
          end if
        end do
      elseif (tp == check_fe(.true.,   4, 4,  4, 0)) then
        !************************************* Quadrangle, Lagrange P1 **************************
        !check 4-node quadrilateral jacobian (see http://mms2.ensmp.fr/ef_paris/technologie/transparents/e_Pathology.pdf)
        do k = 1, elg%nel
          s = dot_product(z(:,elg%mm(3,k))-z(:,elg%mm(1,k)), z(:,elg%mm(2,k))-z(:,elg%mm(1,k))) !projection of a3-a1 in <{a2-a1}>
          t = dot_product(z(:,elg%mm(3,k))-z(:,elg%mm(1,k)), z(:,elg%mm(4,k))-z(:,elg%mm(1,k))) !projection of a3-a1 in <{a4-a1}>
          st = [s, t]
          QJac = [QJ_pos(zz, oz, st, zo, -1, -1), &
                  QJ_pos(zz, oz, st, zo,  1, -1), &
                  QJ_pos(zz, oz, st, zo,  1,  1), &
                  QJ_pos(zz, oz, st, zo, -1,  1)]
          call sfind(QJac, .false., pos) !check whether vertices in (s,t)-plane are well oriented
          if     (size(pos,1) == 2) then
            call swap(elg%mm(pos(1),k), elg%mm(pos(2),k))
          elseif (size(pos,1) == 4) then !in theory, never achieved; it remains for security reasons
            call swap(elg%mm(1,k), elg%mm(2,k))
            call swap(elg%mm(3,k), elg%mm(4,k))
          end if
          !counterclockwise orientation is not ensured by a positive jacobian in the (s,t)-plane: instead, check if 3rd component
          !of the normal vector is positive
          if ((z(1,elg%mm(2,k))-z(1,elg%mm(1,k)))*(z(2,elg%mm(4,k))-z(2,elg%mm(1,k))) - &
              (z(2,elg%mm(2,k))-z(2,elg%mm(1,k)))*(z(1,elg%mm(4,k))-z(1,elg%mm(1,k))) < 0._real64) then
            call swap(elg%mm(1,k), elg%mm(2,k))
            call swap(elg%mm(3,k), elg%mm(4,k))
          end if
        end do
      elseif (tp == check_fe(.true.,   4, 4,  6, 4)) then
        !************************************* Tetrahedron, Lagrange P1 *************************
        do k = 1, elg%nel
          if (.not. detDFT_pos(z(:,elg%mm(1,k)), z(:,elg%mm(2,k)), z(:,elg%mm(3,k)), z(:,elg%mm(4,k)))) &
          call swap(elg%mm(2,k), elg%mm(3,k))
        end do
      elseif (tp == check_fe(.false., 10, 4,  6, 4)) then
        !************************************* Tetrahedron, Lagrange P2 *************************
        do k = 1, elg%nel
          if (.not. detDFT_pos(z(:,elg%mm(1,k)), z(:,elg%mm(2,k)), z(:,elg%mm(3,k)), z(:,elg%mm(4,k)))) &
            call swap(elg%mm(2,k), elg%mm( 3,k))
          if (.not. detDFT_pos(z(:,elg%nn(1,k)), z(:,elg%nn(2,k)), z(:,elg%nn(3,k)), z(:,elg%nn(4,k)))) then
            call swap(elg%nn(2,k), elg%nn( 3,k))
            call swap(elg%nn(5,k), elg%nn( 7,k))
            call swap(elg%nn(9,k), elg%nn(10,k))
          end if
        end do
      elseif (tp == check_fe(.true.,   8, 8, 12, 6)) then
        !************************************* Hexahedron, Lagrange P1 **************************
        !sufficient condition to ensure positive Jacobian for an hexahedron
        !(see http://www.math.udel.edu/~szhang/research/p/subtettest.pdf)
        do k = 1, elg%nel
          check = .true.
          do i = 1, size(Pc, 1)
          if (.not. detDFT_pos(z(:,elg%mm(Pc(i,1),k)), z(:,elg%mm(Pc(i,2),k)), z(:,elg%mm(Pc(i,3),k)), z(:,elg%mm(Pc(i,4),k)))) &
            check = .false.
          end do
          if(.not. check) then
          call info('(module_pmh/reorder_nodes) hexahedron '//trim(string(k))//' does not fulfill sufficient condition to '//&
          &' ensure positive Jacobian: piece '//trim(string(ip))//', group '//trim(string(ig))//'; node order remains unchanged')
          endif
        end do
      else
        call info('(module_pmh/reorder_nodes) reordering of element type '//trim(string(FEDB(pmh%pc(ip)%el(ig)%type)%desc))//&
        &' is not implemented: piece '//trim(string(ip))//', group '//trim(string(ig))//'; node order remains unchanged')
      end if
    end associate
  end do
end do

end subroutine

!-----------------------------------------------------------------------
! swap: swap the value of two variables
!-----------------------------------------------------------------------
subroutine swap(a,b)
integer :: a,b,c
c = a; a = b; b = c
end subroutine

!-----------------------------------------------------------------------
! QJ_pos: check whether the 4-node quadrilateral jacobian is positive
! Only (x,y) coordinates are used
!-----------------------------------------------------------------------
function QJ_pos(a1, a2, a3, a4, r, s) result(res)
real(real64), intent(in) :: a1(2), a2(2), a3(2), a4(2)
integer :: r, s
logical :: res

res = (a4(2)-a2(2))*(a3(1)-a1(1))-(a3(2)-a1(2))*(a4(1)-a2(1))+ &
     ((a3(2)-a4(2))*(a2(1)-a1(1))-(a2(2)-a1(2))*(a3(1)-a4(1)))*r+ &
     ((a4(2)-a1(2))*(a3(1)-a2(1))-(a3(2)-a2(2))*(a4(1)-a1(1)))*s > 0
end function


!-----------------------------------------------------------------------
! detDFT_pos: check whether the tetrahedron jacobian is positive
!-----------------------------------------------------------------------
function detDFT_pos(a1, a2, a3, a4) result(res)
real(real64), intent(in) :: a1(3), a2(3), a3(3), a4(3)
logical :: res

res = (a2(1)-a1(1))*(a3(2)-a1(2))*(a4(3)-a1(3)) + (a2(3)-a1(3))*(a3(1)-a1(1))*(a4(2)-a1(2)) &
    + (a2(2)-a1(2))*(a3(3)-a1(3))*(a4(1)-a1(1)) - (a2(3)-a1(3))*(a3(2)-a1(2))*(a4(1)-a1(1)) &
    - (a2(2)-a1(2))*(a3(1)-a1(1))*(a4(3)-a1(3)) - (a2(1)-a1(1))*(a3(3)-a1(3))*(a4(2)-a1(2)) > 0
end function

!--------------------------------------------------------------------
! cell2node: calculate a node fields from a cell fields
!--------------------------------------------------------------------
subroutine cell2node_real_pmh(pmh)
  type(pmh_mesh), intent(inout) :: pmh
  type(field), allocatable      :: tempfields(:)
  integer, allocatable          :: cell4node(:,:)
  integer :: n_nod_fi,n_cell_fi, maxtopdim
  integer :: ip,ig,ifi,np, nshots
  integer :: ncomp, lnn, nel, i, j, k

  np = 1 ! Only in the first parammeter
  call info('Converting every elementwise field as a nodewise ...')
  ! walk over all pieces
  do ip=1,size(pmh%pc,1)
    associate(pc => pmh%pc(ip))
      n_nod_fi = get_piece_num_fields(pc, 'node')
      n_cell_fi = get_piece_num_fields(pc, 'cell')
      maxtopdim = get_piece_max_top_dim(pc)

      if(n_cell_fi == 0) cycle
      ! Memory allocation for fields
      if(allocated(cell4node)) deallocate(cell4node)
      allocate(cell4node(n_cell_fi, pc%nnod))
      if(n_nod_fi == 0 .and. .not. allocated(pc%fi)) then
        allocate(pc%fi(n_cell_fi))
      else
        if(allocated(tempfields)) then
          do i=1,size(tempfields,1)
            if(allocated(tempfields(i)%val)) deallocate(tempfields(i)%val)
          enddo
          deallocate(tempfields)
        endif
        allocate(tempfields(n_nod_fi+n_cell_fi))
        tempfields(1:n_nod_fi) = pc%fi(1:n_nod_fi)
        call move_alloc(from=tempfields,to=pc%fi)
      endif

      cell4node = 0
      ! walk over all element groups
      do ig=1,size(pc%el,1)
        associate(elg => pc%el(ig))
          if(FEDB(elg%type)%tdim < maxtopdim) cycle
          if(.not. allocated(elg%fi) .or. size(elg%fi,1) /= n_cell_fi) &
            & call error('Wrong number of fields in piece '//trim(string(ip))//&
              & ' group '//trim(string(ig)))
          nel = elg%nel
          lnn = FEDB(elg%type)%lnn

          do ifi=1, n_cell_fi
            if(.not. allocated(elg%fi(ifi)%val)) &
              & call error('Not allocated field '//trim(string(ifi))// &
                & ' in piece '//trim(string(ip))//' group '//trim(string(ig)))
            ncomp = size(elg%fi(ifi)%val,1)
            nshots = size(elg%fi(ifi)%param,1)
            ! Allocate field values, parameters and assign field name.
            if(.not. allocated(pc%fi(n_nod_fi+ifi)%val)) then
              allocate(pc%fi(n_nod_fi+ifi)%val(ncomp,pc%nnod,nshots))
              pc%fi(n_nod_fi+ifi)%val = 0._real64
              pc%fi(n_nod_fi+ifi)%name = trim(elg%fi(ifi)%name)
              call info('  Field: '//trim(elg%fi(ifi)%name))
              if(.not. allocated(pc%fi(n_nod_fi+ifi)%param)) allocate(pc%fi(n_nod_fi+ifi)%param(nshots))
              pc%fi(n_nod_fi+ifi)%param(1:nshots) = elg%fi(ifi)%param(1:nshots)
            endif
            ! Calculate values at nodes
            do k = 1, nel
              do j = 1, lnn
                if(FEDB(elg%type)%nver_eq_nnod) then
                  i = elg%mm(j,k)
                else
                  i = elg%nn(j,k)
                endif
                cell4node(ifi,i) = cell4node(ifi,i) + 1
                do np = 1, nshots
                  pc%fi(n_nod_fi+ifi)%val(:,i,np) = &
                    & pc%fi(n_nod_fi+ifi)%val(:,i,np) + elg%fi(ifi)%val(:,k,np)
                enddo
              enddo
            enddo
            do j = 1, pc%nnod
              do np = 1, nshots
                pc%fi(n_nod_fi+ifi)%val(:,j,np) = pc%fi(n_nod_fi+ifi)%val(:,j,np)/cell4node(ifi,j)
              enddo
            enddo
          enddo
          ! Deallocate cell fields
          if(allocated(elg%fi)) then
            do i=1,size(elg%fi,1)
              if(allocated(elg%fi(i)%val)) deallocate(elg%fi(i)%val)
              if(allocated(elg%fi(i)%param)) deallocate(elg%fi(i)%param)
            enddo
            deallocate(elg%fi)
          endif
        end associate
      enddo
    end associate
  enddo

end subroutine


!--------------------------------------------------------------------
! get_piece_num_fields: returns the number of fields at nodes or on cells or all
!--------------------------------------------------------------------
function get_piece_num_fields(pc, location) result(res)
  type(piece),       intent(in) :: pc
  character(len=*), optional    :: location
  integer                       :: res
  integer                       :: fi_loc !1:node, 2:cell, 3:all
  integer                       :: aux, i

  res = 0
  aux = 0
  if(present(location)) then
    if(trim(adjustl(location)) == 'node') then
      fi_loc = 1 !node
    elseif(trim(adjustl(location)) == 'cell') then
      fi_loc = 2 !cell
    else
      fi_loc = 3 !all
    endif
  else
    fi_loc = 3   !all
  endif

  ! Count node fields
  if(fi_loc==1 .or. fi_loc==3) then
    if(allocated(pc%fi)) then
      res = res+size(pc%fi,1)
    endif
  endif

  ! Count cell fields
  if(fi_loc==2 .or. fi_loc==3) then
    if(allocated(pc%el)) then
      do i=1,size(pc%el,1)
        if(allocated(pc%el(i)%fi)) aux = max(size(pc%el(i)%fi,1),aux)
      enddo
    endif
    res = res + aux
  endif

end function

!--------------------------------------------------------------------
! get_piece_max_top_dim: returns the maximum topological dimension
!--------------------------------------------------------------------
function get_piece_max_top_dim(pc) result(res)
  type(piece),    intent(inout) :: pc
  integer                       :: res
  integer                       :: i

  res = 0

  if(allocated(pc%el)) then
    do i=1, size(pc%el,1)
      res = max(FEDB(pc%el(i)%type)%tdim,res)
    enddo
  endif

end function

!--------------------------------------------------------------------
! get_field_num_shots: returns the number of shots for a given fieldname
!--------------------------------------------------------------------
function get_field_num_shots(pmh,fieldname) result(res)
  type(pmh_mesh),    intent(in) :: pmh
  character(len=*),  intent(in) :: fieldname
  integer                       :: res
  integer                       :: ip, ig, ifi

  res = 0

  do ip=1,size(pmh%pc)
    ! node fields
    if(allocated(pmh%pc(ip)%fi)) then
      do ifi=1, size(pmh%pc(ip)%fi,1)
        if(trim(adjustl(pmh%pc(ip)%fi(ifi)%name)) == trim(adjustl(fieldname))) then
          res = size(pmh%pc(ip)%fi(ifi)%param,1)
          return
        endif
      enddo
    endif

    ! cell fields
    if(allocated(pmh%pc(ip)%el)) then
      do ig=1,size(pmh%pc(ip)%el,1)
        if(allocated(pmh%pc(ip)%el(ig)%fi)) then
          do ifi=1, size(pmh%pc(ip)%el(ig)%fi,1)
            if(trim(adjustl(pmh%pc(ip)%el(ig)%fi(ifi)%name)) == trim(adjustl(fieldname))) then
              res = size(pmh%pc(ip)%el(ig)%fi(ifi)%param,1)
              return
            endif
          enddo
        endif
      enddo
    endif

  enddo

end function

!--------------------------------------------------------------------
! get_num_shots: returns an array with the number of shots of all fields
!--------------------------------------------------------------------
subroutine get_num_shots(pmh, res)
  type(pmh_mesh),    intent(in)    :: pmh
  integer, allocatable, intent(out):: res(:)
  integer                          :: ip, ig, ifi, nfield

  nfield = 0
  do ip=1,size(pmh%pc)
    ! node fields
    if(allocated(pmh%pc(ip)%fi)) then
      do ifi=1, size(pmh%pc(ip)%fi,1)
        nfield = nfield + 1
        call set(res, size(pmh%pc(ip)%fi(ifi)%param,1), nfield, fit=.true.)
      enddo
    endif

    ! cell fields
    if(allocated(pmh%pc(ip)%el)) then
      do ig=1,size(pmh%pc(ip)%el,1)
        if(allocated(pmh%pc(ip)%el(ig)%fi)) then
          do ifi=1, size(pmh%pc(ip)%el(ig)%fi,1)
            nfield = nfield + 1
            call set(res, size(pmh%pc(ip)%el(ig)%fi(ifi)%param,1), nfield, fit=.true.)
          enddo
        endif
      enddo
    endif

  enddo

  if(.not. allocated(res)) allocate(res(0))

end subroutine


!--------------------------------------------------------------------
! remove_coordinate: reduces the space dimension of the mesh removing the chosen coordinate
!--------------------------------------------------------------------
subroutine remove_coordinate(pmh, dim)
  type(pmh_mesh), intent(inout) :: pmh
  integer,        intent(in)    :: dim
  integer                       :: i, j, newdim
  real(real64), allocatable     :: tempz(:,:)

  if(.not. allocated(pmh%pc)) call error('Not allocated mesh')
  call info('Removing component '//string(dim))
  do i=1,size(pmh%pc,1)
    if(pmh%pc(i)%dim<dim) call error('Not enought components. Cannot downgrade space dimension.')
    if(allocated(tempz)) deallocate(tempz)
    newdim = pmh%pc(i)%dim-1
    allocate(tempz(newdim, pmh%pc(i)%nver))
    !tempz(1:newdim,:) = pmh%pc(i)%z([1:dim-1,dim+1:pmh%pc(i)%dim],:)
    tempz(1:newdim,:) = pmh%pc(i)%z((/(j,j=1,dim-1),(j,j=dim+1,pmh%pc(i)%dim)/),:)
    call move_alloc(from=tempz,to=pmh%pc(i)%z)
    pmh%pc(i)%dim = newdim
  enddo
end subroutine


!--------------------------------------------------------------------
! change_pmh_references: changes pmh references
!--------------------------------------------------------------------
subroutine change_pmh_references(pmh, chref)
type(pmh_mesh), intent(inout) :: pmh
integer,        intent(in)    :: chref(:)
integer :: k, ip, ig

call info('Changing references for topological dimension '//trim(string(chref(1)))//'.')
if (allocated(pmh%pc)) then
  do ip = 1, size(pmh%pc,1)
    if (allocated(pmh%pc(ip)%el)) then
      do ig = 1, size(pmh%pc(ip)%el,1)
        associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
          if (FEDB(tp)%tdim == chref(1)) then !same topological dimension
            do k = 2, size(chref,1), 2
              where (elg%ref == chref(k)) elg%ref = chref(k+1)
            end do
          end if
        end associate
      end do
    end if
  end do
end if
end subroutine

end module
