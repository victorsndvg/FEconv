module module_freefem_fcnv
!-----------------------------------------------------------------------
! Module to manage FreeFem++ MSH and MESH meshes
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 20/05/2014
!
! PUBLIC PROCEDURES:
! load_freefem_msh:  load a MSH  file into a PMH structure
! load_freefem_mesh: load a MESH file into a PMH structure
! save_freefem_msh:  save a PMH structure into a MSH  (FreeFem++) file
! save_freefem_mesh: save a PMH structure into a MESH (FreeFem++) file
!-----------------------------------------------------------------------
use basicmod, only: iostat_end, maxpath, error, info, string, int, word_count, lcase, adjustlt, &
                    alloc, dealloc, is_arg, get_post_arg, feed, empty
use module_fe_database_pmh_fcnv, only: FEDB, check_fe
use module_pmh_fcnv, only: pmh_mesh, build_vertices
implicit none

!Private procedures
private :: search_mark

contains

!-----------------------------------------------------------------------
! load_freefem_msh: load a MSH file into a PMH structure
!-----------------------------------------------------------------------
subroutine load_freefem_msh(filename, iu, pmh)
character(*), intent(in) :: filename
integer,      intent(in) :: iu
type(pmh_mesh), intent(inout) :: pmh
integer :: res, i, j, k, ios
character(maxpath) :: cad, str

!allocation
if (allocated(pmh%pc)) then
  deallocate(pmh%pc, stat = res, errmsg = cad)
  if (res /= 0) call error('(module_freefem/load_freefem_msh) Unable to deallocate piece: '//trim(cad))
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_freefem/load_freefem_msh) Unable to allocate piece: '//trim(cad))
else
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_freefem/load_freefem_msh) Unable to allocate piece: '//trim(cad))
end if
allocate(pmh%pc(1)%el(3), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/load_freefem_msh) Unable to allocate groups: '//trim(cad))

!open file
open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('load/open, #'//trim(string(ios)))

associate (m => pmh%pc(1)) !m: current mesh
  !read nver, nel, nel(boundary)
  read (unit=iu, fmt=*, iostat=ios) m%nver, m%el(1)%nel, m%el(2)%nel
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  !determine dim (and therefore, max_tdim) reading the first vertex coordinates
  read (unit=iu, fmt='(a)', iostat=ios) str
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  if (word_count(str) == 3+1) then
    m%dim = 3
    m%el(1)%type = check_fe(.true., 4, 4, 6, 4) !tetrahedra
    m%el(2)%type = check_fe(.true., 3, 3, 3, 0) !triangles
  elseif (word_count(str) == 2+1) then
    m%dim = 2
    m%el(1)%type = check_fe(.true., 3, 3, 3, 0) !triangles
    m%el(2)%type = check_fe(.true., 2, 2, 1, 0) !edges
  else
    call error('(module_freefem/load_freefem_msh) such dimension cannot be saved in a FreeFem++ mesh: '//&
    &trim(string(word_count(str)-1)))
  end if
  !piece
  m%nnod = m%nver
  call alloc(m%z, m%dim, m%nver)
  !element group
  call alloc(m%el(1)%mm, FEDB(m%el(1)%type)%lnv, m%el(1)%nel)
  call alloc(m%el(1)%ref,                        m%el(1)%nel)
  !boundary group
  call alloc(m%el(2)%mm, FEDB(m%el(2)%type)%lnv, m%el(2)%nel)
  call alloc(m%el(2)%ref,                        m%el(2)%nel)
  !node group
  m%el(3)%nel  = m%nver
  m%el(3)%type = check_fe(.true., 1, 1, 0, 0)
  call alloc(m%el(3)%mm, 1, m%el(3)%nel)
  call alloc(m%el(3)%ref,   m%el(3)%nel)
  do j = 1, m%nver; m%el(3)%mm(1,j) = j; end do
  !read z
  backspace(iu)
  if (m%nver > 0) then
    read (unit=iu, fmt=*, iostat=ios) ((m%z(i,j),  i=1,m%dim), m%el(3)%ref(j), j=1,m%nver)
    if (ios /= 0) call error('(module_freefem/load_freefem_msh) unable to read z: #'//trim(string(ios)))
  end if  
  !read el(1)%mm
  if (m%el(1)%nel > 0) then
    read (unit=iu, fmt=*, iostat=ios) ((m%el(1)%mm(i,k),  i=1,FEDB(m%el(1)%type)%lnv), m%el(1)%ref(k), k=1,m%el(1)%nel)
    if (ios /= 0) call error('(module_freefem/load_freefem_msh) unable to read el(1)%mm: #'//trim(string(ios)))
  end if
  !read el(2)%mm
  if (m%el(2)%nel > 0) then
    read (unit=iu, fmt=*, iostat=ios) ((m%el(2)%mm(i,k),  i=1,FEDB(m%el(2)%type)%lnv), m%el(2)%ref(k), k=1,m%el(2)%nel)
    if (ios /= 0) call error('(module_freefem/load_freefem_msh) unable to read el(2)%mm: #'//trim(string(ios)))
  end if
end associate
call build_vertices(pmh)
end subroutine

!-----------------------------------------------------------------------
! load_freefem_mesh: load a MESH file into a PMH structure
!-----------------------------------------------------------------------
subroutine load_freefem_mesh(filename, iu, pmh)
character(*), intent(in) :: filename
integer,      intent(in) :: iu
type(pmh_mesh), intent(inout) :: pmh
integer :: res, i, j, k, ios
character(maxpath) :: cad

!allocation
if (allocated(pmh%pc)) then
  deallocate(pmh%pc, stat = res, errmsg = cad)
  if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to deallocate piece: '//trim(cad))
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to allocate piece: '//trim(cad))
else
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to allocate piece: '//trim(cad))
end if
allocate(pmh%pc(1)%el(3), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to allocate groups: '//trim(cad))
!open file
open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('load/open, #'//trim(string(ios)))
!check if there are tetrahedra
res = search_mark(iu, 'Tetrahedra')
if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to find mark Tetrahedra: #'//trim(string(res)))
rewind(iu)

associate (m => pmh%pc(1)) !m: current mesh
  m%dim = 3
  m%el(1)%type = check_fe(.true., 4, 4, 6, 4) !tetrahedra
  m%el(2)%type = check_fe(.true., 3, 3, 3, 0) !triangles
  !read nver and z
  res = search_mark(iu, 'Vertices')
  if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to find mark Vertices: #'//trim(string(res)))
  read (unit=iu, fmt=*, iostat=ios) m%nver
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  m%nnod = m%nver
  call alloc(m%z, m%dim, m%nver)
  m%el(3)%nel  = m%nver !node group
  m%el(3)%type = check_fe(.true., 1, 1, 0, 0)
  call alloc(m%el(3)%mm, 1, m%el(3)%nel)
  call alloc(m%el(3)%ref,   m%el(3)%nel)
  do j = 1, m%nver; m%el(3)%mm(1,j) = j; end do
  read (unit=iu, fmt=*, iostat=ios) ((m%z(i,j),  i=1,m%dim), m%el(3)%ref(j), j=1,m%nver)
  if (ios /= 0) call error('(module_freefem/load_freefem_mesh) unable to read z: #'//trim(string(ios)))
  !read el(1)
  rewind(iu)
  res = search_mark(iu, 'Tetrahedra')
  if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to find mark Tetrahedra: #'//trim(string(res)))
  read (unit=iu, fmt=*, iostat=ios) m%el(1)%nel
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  call alloc(m%el(1)%mm, FEDB(m%el(1)%type)%lnv, m%el(1)%nel)
  call alloc(m%el(1)%ref,                        m%el(1)%nel)
  read (unit=iu, fmt=*, iostat=ios) ((m%el(1)%mm(i,k),  i=1,FEDB(m%el(1)%type)%lnv), m%el(1)%ref(k), k=1,m%el(1)%nel)
  if (ios /= 0) call error('(module_freefem/load_freefem_mesh) unable to read el(1): #'//trim(string(ios)))
  !read el(2)
  rewind(iu)
  res = search_mark(iu, 'Triangles')
  if (res /= 0) call error('(module_freefem/load_freefem_mesh) Unable to find mark Triangles: #'//trim(string(res)))
  read (unit=iu, fmt=*, iostat=ios) m%el(2)%nel
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  call alloc(m%el(2)%mm, FEDB(m%el(2)%type)%lnv, m%el(2)%nel)
  call alloc(m%el(2)%ref,                        m%el(2)%nel)
  read (unit=iu, fmt=*, iostat=ios) ((m%el(2)%mm(i,k),  i=1,FEDB(m%el(2)%type)%lnv), m%el(2)%ref(k), k=1,m%el(2)%nel)
  if (ios /= 0) call error('(module_freefem/load_freefem_mesh) unable to read el(2): #'//trim(string(ios)))
end associate
call build_vertices(pmh)
end subroutine

!-----------------------------------------------------------------------
! save_freefem: save a PMH structure into a MSH (FreeFem++) file
!
! pmh is deallocated while variables are being saved
!-----------------------------------------------------------------------
subroutine save_freefem_msh(outfile, iu, pmh)
character(*),   intent(in)    :: outfile
integer,        intent(in)    :: iu
type(pmh_mesh), intent(inout) :: pmh

integer :: i, ipp, ip, ig, k, j, type_by_tdim(0:3), prev_max_tdim, res, max_tdim, ios, nel, ntri, nver, dim, el_type, bnd_type
integer, allocatable :: piece2save(:), el_piece(:), bnd_piece(:), nver_piece(:)
character(maxpath) :: str, cad

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
  call info('(module_freefem/save_freefem) option -glue not implemented yet')
end if

!testing and calculation of max_tdim
type_by_tdim  = 0 !store the type of element for each topological dimension
prev_max_tdim = 0 !store the maximal topological dimension
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_freefem/save_freefem) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (tp /= check_fe(.true., 4, 4, 6, 4) .and. tp /= check_fe(.true., 3, 3, 3, 0) .and. tp /= check_fe(.true., 2, 2, 1, 0)) then
        call info('(module_freefem/save_freefem) element type '//trim(FEDB(tp)%desc)//' found; those elements cannot be saved'//&
        &' in FreeFem++ format and they will be discarded')
        cycle
      end if
      !check whether there is only one type of element for each topological dimension
      if (type_by_tdim( FEDB(tp)%tdim ) == 0) then
        type_by_tdim( FEDB(tp)%tdim ) = tp
      elseif (type_by_tdim( FEDB(tp)%tdim ) /= tp) then
        call error('(module_freefem/save_freefem) more that one type of element is defined for the same topological dimension: '//&
        &string(type_by_tdim(FEDB(tp)%tdim))//', '//string(tp)//'; unable to convert to FreeFem++')
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
    call error('(module_freefem/save_freefem) there are pieces with different maximal topological dimension; unable to convert '//&
    &'to FreeFem++')
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
    call error('(module_freefem/save_freefem) different coordinates dimensions in different pieces: '//&
    &string(pmh%pc(ip)%dim)//', '//string(dim)//'; unable to convert to FreeFem++')
  end if
end do
if (max_tdim == 3) then !tetrahedra
  el_type  = check_fe(.true., 4, 4, 6, 4)
  bnd_type = check_fe(.true., 3, 3, 3, 0)
elseif (max_tdim == 2) then !triangles
  el_type  = check_fe(.true., 3, 3, 3, 0)
  bnd_type = check_fe(.true., 2, 2, 1, 0)
else
  call error('(module_freefem/save_freefem) such topological dimension cannot be saved in a FreeFem++ mesh: '//&
  &trim(string(max_tdim)))
end if

!store variables nnod, nver, nel for selected pieces
if (allocated(el_piece)) deallocate(el_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(el_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable el_piece: '//trim(cad))
if (allocated(bnd_piece)) deallocate(bnd_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(bnd_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable bnd_piece: '//trim(cad))
if (allocated(nver_piece)) deallocate(nver_piece); !nver_piece(ipp): global numbering for the last vertex of piece #ipp
allocate(nver_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable nver_piece: '//trim(cad))
el_piece(0) = 0; bnd_piece(0) = 0; nver_piece(0) = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  nver_piece(ipp) = nver_piece(ipp-1) + pmh%pc(ipp)%nver
  el_piece(ipp)   =   el_piece(ipp-1)
  bnd_piece(ipp)  =  bnd_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == el_type) then
         el_piece(ipp)  =   el_piece(ipp) + elg%nel
      elseif (elg%type == bnd_type) then
         bnd_piece(ipp) =  bnd_piece(ipp) + elg%nel
      end if
    end associate
  end do
end do

nel  =   el_piece(size(piece2save,1))
ntri =  bnd_piece(size(piece2save,1))
nver = nver_piece(size(piece2save,1))
!store variables nver, nel, ntri for selected pieces
open (unit=iu, file=outfile, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('save/open, #'//trim(string(ios)))
call feed(iu, string(nver)); call feed(iu, string(nel)); call feed(iu, string(ntri));  call empty(iu)
!store vertex coordinates and reference 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do j = 1, pmh%pc(ip)%nver
    if (pmh%pc(ip)%dim == max_tdim) then
      do i = 1, pmh%pc(ip)%dim;   call feed(iu, string(pmh%pc(ip)%z(i,j))); end do
    elseif (max_tdim == 3 .and. pmh%pc(ip)%dim < max_tdim) then
      do i = 1, pmh%pc(ip)%dim;   call feed(iu, string(pmh%pc(ip)%z(i,j))); end do
      do i = 1, 3-pmh%pc(ip)%dim; call feed(iu, string(0.));                end do
    elseif (max_tdim == 2 .and. pmh%pc(ip)%dim > max_tdim) then
      do i = 1, max_tdim;         call feed(iu, string(pmh%pc(ip)%z(i,j))); end do
    end if
    call feed(iu, string(0))
    call empty(iu)
  end do
end do
!store element conectivities and references
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == el_type) then
        do k = 1, elg%nel
          do i = 1, FEDB(el_type)%lnv; call feed(iu, string(nver_piece(ipp-1)+elg%mm(i,k))); end do
          call feed(iu, string(elg%ref(k)))
          call empty(iu)
        end do
      end if
    end associate
  end do
end do
!store boundary conectivities and references
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == bnd_type) then
        do k = 1, elg%nel
          do i = 1, FEDB(bnd_type)%lnv; call feed(iu, string(nver_piece(ipp-1)+elg%mm(i,k))); end do
          call feed(iu, string(elg%ref(k)))
          call empty(iu)
        end do
      end if
    end associate
  end do
end do
close(iu)
end subroutine

!-----------------------------------------------------------------------
! save_freefem: save a PMH structure into a MESH (FreeFem++) file
!
! pmh is deallocated while variables are being saved
!-----------------------------------------------------------------------
subroutine save_freefem_mesh(outfile, iu, pmh)
character(*),   intent(in)    :: outfile
integer,        intent(in)    :: iu
type(pmh_mesh), intent(inout) :: pmh

integer :: i, ipp, ip, ig, k, j, type_by_tdim(0:3), prev_max_tdim, res, max_tdim, ios, nel, ntri, nver, dim
integer, allocatable :: piece2save(:), el_piece(:), bnd_piece(:), nver_piece(:)
character(maxpath) :: str, cad

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
  call info('(module_freefem/save_freefem) option -glue not implemented yet')
end if

!testing and calculation of max_tdim
type_by_tdim  = 0 !store the type of element for each topological dimension
prev_max_tdim = 0 !store the maximal topological dimension
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_freefem/save_freefem) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (tp /= check_fe(.true., 4, 4, 6, 4) .and. tp /= check_fe(.true., 3, 3, 3, 0)) then
        call info('(module_freefem/save_freefem) element type '//trim(FEDB(tp)%desc)//' found; those elements cannot be saved'//&
        &' in FreeFem++ format and they will be discarded')
        cycle
      end if
      !check whether there is only one type of element for each topological dimension
      if (type_by_tdim( FEDB(tp)%tdim ) == 0) then
        type_by_tdim( FEDB(tp)%tdim ) = tp
      elseif (type_by_tdim( FEDB(tp)%tdim ) /= tp) then
        call error('(module_freefem/save_freefem) more that one type of element is defined for the same topological dimension: '//&
        &string(type_by_tdim(FEDB(tp)%tdim))//', '//string(tp)//'; unable to convert to FreeFem++')
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
    call error('(module_freefem/save_freefem) there are pieces with different maximal topological dimension; unable to convert '//&
    &'to FreeFem++')
  end if
end do
if (max_tdim < 3) call error('(module_freefem/save_freefem) only the conversion to a thetrahedral FreeFem++ mesh is implemented.')
!testing and calculation of dim
dim = -1
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  !check whether there is only one coordinates dimension for all pieces
  if (dim == -1) then
    dim = pmh%pc(ip)%dim
  elseif (pmh%pc(ip)%dim /= dim) then
    call error('(module_freefem/save_freefem) different coordinates dimensions in different pieces: '//&
    &string(pmh%pc(ip)%dim)//', '//string(dim)//'; unable to convert to FreeFem++')
  end if
end do

!store variables nnod, nver, nel for selected pieces
if (allocated(el_piece)) deallocate(el_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(el_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable el_piece: '//trim(cad))
if (allocated(bnd_piece)) deallocate(bnd_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(bnd_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable bnd_piece: '//trim(cad))
if (allocated(nver_piece)) deallocate(nver_piece); !nver_piece(ipp): global numbering for the last vertex of piece #ipp
allocate(nver_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable nver_piece: '//trim(cad))
el_piece(0) = 0; bnd_piece(0) = 0; nver_piece(0) = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  nver_piece(ipp) = nver_piece(ipp-1) + pmh%pc(ipp)%nver
  el_piece(ipp)  =  el_piece(ipp-1)
  bnd_piece(ipp)  =  bnd_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == check_fe(.true., 4, 4, 6, 4)) then
         el_piece(ipp) =  el_piece(ipp) + elg%nel
      elseif (elg%type == check_fe(.true., 3, 3, 3, 0)) then
         bnd_piece(ipp) =  bnd_piece(ipp) + elg%nel
      end if
    end associate
  end do
end do
nel  =  el_piece(size(piece2save,1))
ntri =  bnd_piece(size(piece2save,1))
nver = nver_piece(size(piece2save,1))

!store variables nver, nel, ntri for selected pieces
open (unit=iu, file=outfile, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('save/open, #'//trim(string(ios)))
!call feed(iu, string(nver)); call feed(iu, string(nel)); call feed(iu, string(ntri));  call empty(iu)

write(iu, '(a/)') 'MeshVersionFormatted 1'
write(iu, '(a)')  'Vertices'
write(iu, '(a)')  trim(string(nver))

!store vertex coordinates and reference 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do j = 1, pmh%pc(ip)%nver
    do i = 1, pmh%pc(ip)%dim; call feed(iu, string(pmh%pc(ip)%z(i,j))); end do
    do i = 1, 3-dim;          call feed(iu, string(0.));                end do
    call feed(iu, string(0))
    call empty(iu)
  end do
end do

!store tet conectivities and references
write(iu, '(/a)') 'Tetrahedra'
write(iu, '(a)')  trim(string(nel))
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == check_fe(.true., 4, 4, 6, 4)) then
        do k = 1, elg%nel
          do i = 1, 4; call feed(iu, string(nver_piece(ipp-1)+elg%mm(i,k))); end do
          call feed(iu, string(elg%ref(k)))
          call empty(iu)
        end do
      end if
    end associate
  end do
end do

!store tri conectivities and references
write(iu, '(/a)') 'Triangles'
write(iu, '(a)')  trim(string(ntri))
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == check_fe(.true., 3, 3, 3, 0)) then
        do k = 1, elg%nel
          do i = 1, 3; call feed(iu, string(nver_piece(ipp-1)+elg%mm(i,k))); end do
          call feed(iu, string(elg%ref(k)))
          call empty(iu)
        end do
      end if
    end associate
  end do
end do

write(iu, '(/a)') 'End'
close(iu)
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! search_mark: searches for a mark
! RETURN: 0 if the mark is found; iostat_end if end-of-file is found;
! non-zero otherwise
!-----------------------------------------------------------------------
function search_mark(id, mark) result(res)
integer,      intent(in) :: id
character(*), intent(in) :: mark
character(maxpath) :: str
integer :: res

do
  read (id, '(a)', iostat=res) str
  if (res == iostat_end) return !end-of-file found
  if (res /= 0) call error('(module_freefem/search/read) error while reading mark "'//trim(mark)//'", #'//trim(string(res)))
  if (adjustlt(lcase(str)) == lcase(mark)) return
end do
end function

end module
