module module_pmh
!-----------------------------------------------------------------------
! Module to manage piecewise meshes (PMH)
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 19/01/2014
!
! PUBLIC PROCEDURES:
!   pmh2mfm: convert a PMH structure into a MFM one
!   mfm2pmh: convert a MFM mesh into a PMH structure
!   build_vertices: build vertex connectivities and coordinates from node information (P1, P2 only)
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
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_convers, only: string, int
use module_alloc, only: alloc, dealloc, reduce, find_first, find_row_sorted, sort, insert, insert_row_sorted, bsearch
use module_args, only: is_arg, get_post_arg
use module_fe_database_pmh, only : FEDB, check_fe
implicit none

!Types
type elgroup
  integer              :: nel  = 0 !total number of elements
  integer              :: type = 0 !element type (one of those defined in module_eltype)
  integer, allocatable :: nn(:,:)  !global numbering of nodes
  integer, allocatable :: mm(:,:)  !global numbering of vertices
  integer, allocatable :: ref(:)   !reference numbering
end type

type piece
  integer                    :: nnod = 0 !total number of nodes
  integer                    :: nver = 0 !total number of vertices
  integer                    :: dim  = 0 !space dimension of the node/vertex coordinates
  real(real64),  allocatable :: z(:,:)   !vertex coordinates
  type(elgroup), allocatable :: el(:)    !element groups  
end type

type pmh_mesh
  type(piece), allocatable :: pc(:) !pieces that compose the mesh
end type  

contains

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
prev_max_tdim, res, max_tdim, valid_fe(12)
integer, allocatable :: piece2save(:), ref(:,:), tmp_vf(:), nel_piece(:), nnod_piece(:), nver_piece(:)
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
            check_fe(.true.,   8, 8, 12, 6)]   !Hexahedron, Lagrange P1

!check piece(s) to be saved
if (is_arg('-p')) then 
  str = get_post_arg('-p')
  if (str == '0') then !save all pieces
    call alloc(piece2save, size(pmh%pc,1))
    piece2save = [(i, i=1, size(pmh%pc,1))]
  else !save a single piece, indicated in -p
    call alloc(piece2save, 1)
    piece2save(1) = int(str)
  end if
else !save only the first piece
  call alloc(piece2save, 1)
  piece2save(1) = 1
end if
if (is_arg('-glue')) then 
  call error('(module_pmh/pmh2mfm) option -glue not implemented yet')
end if

!testing and calculation of max_tdim
type_by_tdim  = 0 !store the type of element for each topological dimension
prev_max_tdim = 0 !store the maximal topological dimension
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
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
  nnod_piece(ipp) = nnod_piece(ipp-1) + pmh%pc(ipp)%nnod
  nver_piece(ipp) = nver_piece(ipp-1) + pmh%pc(ipp)%nver
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
            call insert_row_sorted(ref, tmp_2d, used=n, fit=[.false.,.true.])
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
          pos = find_first(ref(1,:), mm(j,k))
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
            call insert_row_sorted(ref, tmp_3d, used=n, fit=[.false.,.true.])
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
          pos = find_row_sorted(ref, tmp_2d, n)
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
              call insert_row_sorted(ref, tmp_vf, used=n, fit=[.false.,.true.])
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
            pos = find_row_sorted(ref, tmp_vf, n)
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
  z(1:dim, nver_piece(ipp-1)+1:nver_piece(ipp)) = pmh%pc(ipp)%z
  deallocate(pmh%pc(ipp)%z)
end do
end subroutine

!-----------------------------------------------------------------------
! mfm2pmh: convert a MFM mesh into a PMH structure
!
! mfm is deallocated while PMH variables are being allocated
!-----------------------------------------------------------------------
subroutine mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh)
type(pmh_mesh), intent(inout) :: pmh
integer, intent(inout) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf 
integer, allocatable   :: nn(:,:), mm(:,:), nrv(:,:), nra(:,:), nrc(:,:), nsd(:)
real(real64), allocatable :: z(:,:)
integer :: res, tp, nelg, iel, ic, j, k
character(maxpath) :: cad

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
call move_alloc(from=z, to=pmh%pc(1)%z)

!nelg: calculate the number of element groups to create
tp = check_fe(nver==nnod, lnn, lnv, lne, lnf)
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
  call move_alloc(from=nn,  to=pmh%pc(1)%el(nelg)%nn)
  call move_alloc(from=mm,  to=pmh%pc(1)%el(nelg)%mm)
  call move_alloc(from=nsd, to=pmh%pc(1)%el(nelg)%ref)
end associate
iel = nelg - 1

!group for faces
if (FEDB(tp)%tdim > 2) then
  if (maxval(nrc) > 0) then
    associate (melg => pmh%pc(1)%el(nelg), elg => pmh%pc(1)%el(iel)) !melg: max. tdim group, elg: current group
      elg%nel = count(nrc > 0)
      elg%type = FEDB(melg%type)%f_type
      !nn is not relevant for groups without maximal topological dimension
      call alloc(elg%mm, FEDB(elg%type)%lnv, elg%nel) 
      call alloc(elg%ref, elg%nel) 
      ic = 1
      do k = 1, melg%nel
        do j = 1, FEDB(melg%type)%lnf
          if (nrc(j,k) /= 0) then
            elg%mm(:, ic) = pmh%pc(1)%el(nelg)%mm(FEDB(melg%type)%face(1:FEDB(elg%type)%lnv, j), k)
            elg%ref(ic)   = nrc(j,k)
            ic = ic + 1
          end if  
        end do
      end do
    end associate
    iel = iel - 1
  end if
  call dealloc(nrc)
end if

!group for edges
if (FEDB(tp)%tdim > 1) then
  if (maxval(nra) > 0) then
    associate (melg => pmh%pc(1)%el(nelg), elg => pmh%pc(1)%el(iel)) !melg: max. tdim group, elg: current group
      elg%nel = count(nra > 0)
      elg%type = FEDB(melg%type)%e_type
      !nn is not relevant for groups without maximal topological dimension
      call alloc(elg%mm, FEDB(elg%type)%lnv, elg%nel) 
      call alloc(elg%ref, elg%nel) 
      ic = 1
      do k = 1, melg%nel
        do j = 1, FEDB(melg%type)%lne
          if (nra(j,k) /= 0) then
            elg%mm(:, ic) = pmh%pc(1)%el(nelg)%mm(FEDB(melg%type)%edge(1:FEDB(elg%type)%lnv, j), k)
            elg%ref(ic)   = nra(j,k)
            ic = ic + 1
          end if  
        end do
      end do
    end associate
    iel = iel - 1
  end if
  call dealloc(nra)
end if

!group for vertices
if (FEDB(tp)%tdim > 0) then
  if (maxval(nrv) > 0) then
    associate (melg => pmh%pc(1)%el(nelg), elg => pmh%pc(1)%el(iel)) !melg: max. tdim group, elg: current group
      elg%nel = count(nrv > 0)
      elg%type = FEDB(melg%type)%v_type
      !nn is not relevant for groups without maximal topological dimension
      call alloc(elg%mm, FEDB(elg%type)%lnv, elg%nel) 
      call alloc(elg%ref, elg%nel) 
      ic = 1
      do k = 1, melg%nel
        do j = 1, FEDB(melg%type)%lnv
          if (nrv(j,k) /= 0) then
            elg%mm(:, ic) = pmh%pc(1)%el(nelg)%mm(FEDB(melg%type)%edge(1:FEDB(elg%type)%lnv, j), k)
            elg%ref(ic)   = nrv(j,k)
            ic = ic + 1
          end if  
        end do
      end do
    end associate
    iel = iel - 1
  end if
  call dealloc(nrv)
end if

end subroutine

!-----------------------------------------------------------------------
! build_vertices: build vertex connectivities and coordinates from node information
!
! only valid for Lagrange P1 and P2 elements
! it is assumed that vertices are the first nodes
! it is assumed that global numbering as vertex is less or equal than global numbering as node
!-----------------------------------------------------------------------
subroutine build_vertices(pmh)
type(pmh_mesh), intent(inout) :: pmh
integer, allocatable :: vert2node(:), node2vert(:)
integer :: nv2d, pos, maxv, i, j, k, ig, ip, valid_fe(10)
logical nver_eq_nnod

!valid elements types to save a MFM mesh 
valid_fe = [check_fe(.true.,   1, 1,  0, 0), & !Node
            check_fe(.true.,   2, 2,  1, 0), & !Edge, Lagrange P1 
            check_fe(.false.,  3, 2,  1, 0), & !Edge, Lagrange P2
            check_fe(.true.,   3, 3,  3, 0), & !Triangle, Lagrange P1
            check_fe(.true.,   4, 4,  4, 0), & !Quadrangle, Lagrange P1
            check_fe(.false.,  6, 3,  3, 0), & !Triangle, Lagrange P2
            check_fe(.true.,   4, 4,  6, 4), & !Tetrahedron, Lagrange P1
            check_fe(.false., 10, 4,  6, 4), & !Tetrahedron, Lagrange P2
            check_fe(.true.,   8, 8, 12, 6), & !Hexahedron, Lagrange P1
            check_fe(.true.,   6, 6,  9, 5)]   !Wedge, Lagrange P1

!check if all the elements are Lagrange P1
nver_eq_nnod = .true.
do ip = 1, size(pmh%pc,1)
  do ig = 1, size(pmh%pc(ip)%el,1)
    if (find_first(valid_fe, pmh%pc(ip)%el(ig)%type) == 0) then
      call info('(module_pmh/build_vertices) element type '//trim(FEDB(pmh%pc(ip)%el(ig)%type)%desc)//' found; those elements'//&
      &' cannot be used to construct vertex information and they will be discarded')
      cycle
    end if  
    nver_eq_nnod = nver_eq_nnod .and. FEDB(pmh%pc(ip)%el(ig)%type)%nver_eq_nnod
  end do
end do
!all elements are Lagrange P1: vertex information is the same than node information
if (nver_eq_nnod) then
  do ip = 1, size(pmh%pc,1)
    if (.not.allocated(pmh%pc(ip)%z)) call error('(module_pmh/build_vertices) z is not allocated: piece '//trim(string(ip))//&
    &'; unable to build vertices')
    do ig = 1, size(pmh%pc(ip)%el,1)
      if (pmh%pc(ip)%nver == 0) pmh%pc(ip)%nver = pmh%pc(ip)%nnod
      associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type) !elg: current group, tp: element type
        if (find_first(valid_fe, tp) == 0) cycle
        if (allocated(elg%nn) .and. .not.allocated(elg%mm)) then
          call move_alloc(from=elg%nn, to=elg%mm)
        elseif (.not.allocated(elg%nn) .and. .not.allocated(elg%mm)) then
          call error('(module_pmh/build_vertices) neither mm nor nn are not allocated: piece '//trim(string(ip))//&
          &', group '//trim(string(ig))//'; unable to build vertices')
        end if
      end associate      
    end do
  end do
  return
end if

!transformation for Lagrange P1/P2 elements
do ip = 1, size(pmh%pc,1)
  if (.not.allocated(pmh%pc(ip)%z)) call error('(module_pmh/build_vertices)  z is not allocated: piece '//trim(string(ip))//&
  &'; unable to build vertices')
  nv2d = 0
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type) !elg: current group, tp: element type
      if (find_first(valid_fe, tp) == 0) cycle !proceed only with Lagrange P1 or P2 elements
      if (.not.allocated(elg%nn)) call error('(module_pmh/build_vertices) nn is not allocated: piece '//trim(string(ip))//&
      &', group '//trim(string(ig))//'; unable to build vertices')
      !vert2node: stores global numbering of vertices (numbered as nodes) for all groups in a piece; 
      do k = 1, elg%nel
        do i = 1, FEDB(tp)%lnv
          if (elg%nn(i,k) == 0) call error('(module_pmh/build_vertices) node numbering is zero: piece '//trim(string(ip))//&
          &', group '//trim(string(ig))//', node '//trim(string(i))//', element '//trim(string(k))//'; unable to build vertices')
          pos = bsearch(vert2node, elg%nn(i,k), nv2d) 
          if (pos < 0) then
            call insert(vert2node, elg%nn(i,k), -pos, nv2d, fit=.false.)
            nv2d = nv2d + 1 
          endif
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
  allocate(node2vert(maxv))
  node2vert = 0
  do i = 1, size(vert2node, 1)
    if (vert2node(i) == 0) cycle
    node2vert(vert2node(i)) = i
  enddo
  !vertices renumbering
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type) !elg: current group, tp: element type
      if (find_first(valid_fe, tp) == 0) cycle
      call alloc(elg%mm, FEDB(tp)%lnv, elg%nel)
      do k = 1, elg%nel
        do i = 1, FEDB(tp)%lnv
          elg%mm(i,k) = node2vert(elg%nn(i,k))
        end do
      end do
    end associate
  end do
  !vertices coordinates: assume that i <= vert2node(i); thus z can be overwritten
  associate(pi => pmh%pc(ip))
    do j = 1, pi%nver
      pi%z(1:pi%dim, j) = pi%z(1:pi%dim, vert2node(j))
    end do
    call reduce(pi%z, pi%dim, pi%nver)
  end associate
  call dealloc(vert2node)
  call dealloc(node2vert)
end do
end subroutine

end module

! PASOS A DAR:
!   Definición de pmh
!   Creacion de funciones para procesar las coordenadas, el nn, etc
!   Opcionalmente, creación de funciones para procesar refs, orientación, etc.
!   Conversor a/de mfm
!   Conversor a vtu
!   Uso en comsol
!   hay que hacer subrutinas para dado p2, construir P1 en todos los elementos)
