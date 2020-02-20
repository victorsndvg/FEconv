module module_gmsh_fcnv
!-----------------------------------------------------------------------
! Module to manage Gmsh MSH meshes
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 26/05/2014
!
! PUBLIC PROCEDURES:
! load_gmsh: load a Gmsh ASCII MSH  file into a PMH structure
! save_gmsh: save a Gmsh MSH file from a PMH structure
!-----------------------------------------------------------------------
use basicmod, only: real64, iostat_end, maxpath, error, info, int, string, lcase, adjustlt, is_arg, get_post_arg, &
                    alloc, dealloc, insert_sorted, set, find_sorted, reduce, find_first
use module_fe_database_pmh_fcnv, only: FEDB, check_fe
use module_pmh_fcnv, only: pmh_mesh, build_vertices, build_node_coordinates
implicit none

!Private procedures
private :: search_mark

contains

!-----------------------------------------------------------------------
! load_gmsh: load a Gmsh MSH file into a PMH structure
!-----------------------------------------------------------------------
subroutine load_gmsh(filename, iu, pmh)
character(*), intent(in) :: filename
integer,      intent(in) :: iu
type(pmh_mesh), intent(inout) :: pmh
integer :: res, i, j, k, ios, nel, elt, nelt, id4gmsh(19), xxx, ntags, tag(100), tmp(100)
integer, allocatable :: el_type(:), inv_el_type(:), new4old(:), old4new(:)
character(maxpath) :: cad
real :: version
!From http://gmsh.info//doc/texinfo/gmsh.html#MSH-file-format: elm-type: defines the geometrical type of the n-th element
character(*), parameter :: gmsh_elm_type(19) = [ &
'2-node line                                                                                                                  ', &
'3-node triangle                                                                                                              ', &
'4-node quadrangle                                                                                                            ', &
'4-node tetrahedron                                                                                                           ', &
'8-node hexahedron                                                                                                            ', &
'6-node prism                                                                                                                 ', &
'5-node pyramid                                                                                                               ', &
'3-node second order line (2 nodes associated with the vertices and 1 with the edge)                                          ', &
'6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)                                     ', &
'9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)                  ', &
'10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)                                 ', &
'27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with volume)', &
'18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)        ', &
'14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)       ', &
'1-node point                                                                                                                 ', &
'8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)                                   ', &
'20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)                                 ', &
'15-node second order prism (6 nodes associated with the vertices and 9 with the edges)                                       ', &
'13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)                                     ']

!id4gmsh("gmsh element id") = "pmh element id"
id4gmsh     = 0
id4gmsh(15) = check_fe(.true.,      1, 1,  0, 0) ! 1 Vertex                            
id4gmsh( 1) = check_fe(.true.,      2, 2,  1, 0) ! 2 Edge, Lagrange P1                 
id4gmsh( 8) = check_fe(.false.,     3, 2,  1, 0) ! 3 Edge, Lagrange P2                 
id4gmsh( 2) = check_fe(.true.,      3, 3,  3, 0) ! 4 Triangle, Lagrange P1             
id4gmsh( 9) = check_fe(.false.,     6, 3,  3, 0) ! 5 Triangle, Lagrange P2             
!id4gmsh( ) = check_fe(.false.,     3, 3,  3, 0) ! 6 Triangle, Raviart-Thomas (edge)   
id4gmsh( 3) = check_fe(.true.,      4, 4,  4, 0) ! 7 Quadrangle, Lagrange P1           
id4gmsh(16) = check_fe(.false.,     8, 4,  4, 0) ! 8 Quadrangle, Lagrange P2           
id4gmsh( 4) = check_fe(.true.,      4, 4,  6, 4) ! 9 Tetrahedron, Lagrange P1          
id4gmsh(11) = check_fe(.false.,    10, 4,  6, 4) !10 Tetrahedron, Lagrange P2          
!id4gmsh( ) = check_fe(.false.,     4, 4,  6, 4) !11 Tetrahedron, Raviart-Thomas (face)
!id4gmsh( ) = check_fe(.false.,     6, 4,  6, 4) !12 Tetrahedron, Nedelec (edge)       
id4gmsh( 5) = check_fe(.true.,      8, 8, 12, 6) !13 Tetrahedron, Nedelec 2 (edge)     
id4gmsh(17) = check_fe(.false.,    20, 8, 12, 6) !14 Hexahedron, Lagrange P1           
id4gmsh( 6) = check_fe(.true.,      6, 6,  9, 5) !15 Hexahedron, Lagrange P2           

!allocation
if (allocated(pmh%pc)) then
  deallocate(pmh%pc, stat = res, errmsg = cad)
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to deallocate piece: '//trim(cad))
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to allocate piece: '//trim(cad))
else
  allocate(pmh%pc(1), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to allocate piece: '//trim(cad))
end if
!open file
open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('load/open, #'//trim(string(ios)))
associate (m => pmh%pc(1)) !m: current mesh
  !version
  res = search_mark(iu, '$MeshFormat')
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to find mark $MeshFormat: #'//trim(string(res)))
  read (unit=iu, fmt=*, iostat=ios) version
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  if (version < 2. .or. 3. <= version) call error('(module_gmsh/load_gmsh) Unsupported Gmsh mesh format version '//&
  &trim(string(int(version)))//'. Please, export a "Version 2 ASCII" Gmsh mesh.')
  !nodes
  res = search_mark(iu, '$Nodes')
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to find mark $Nodes: #'//trim(string(res)))
  m%dim = 3
  read (unit=iu, fmt=*, iostat=ios) m%nnod
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  call alloc(m%z, m%dim, m%nnod)
  call alloc(old4new,    m%nnod)
  read (unit=iu, fmt=*, iostat=ios) (old4new(j), (m%z(i,j),  i=1,m%dim), j=1,m%nnod) !j is the new node index, old4new(j), the old
  if (ios /= 0) call error('(module_gmsh/load_gmsh) unable to read z: #'//trim(string(ios)))
  !invert old4new: new4old
  call alloc(new4old, maxval(old4new))
  do j = 1, m%nnod
    new4old(old4new(j)) = j
  end do
  !element types: el_type
  res = search_mark(iu, '$Elements')
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to find mark $Elements: #'//trim(string(res)))
  read (unit=iu, fmt=*, iostat=ios) nel
  if (ios /= 0) call error('read (str), #'//trim(string(ios)))
  call alloc(el_type, nel)
  nelt = 0
  do k = 1, nel
    read (unit=iu, fmt=*, iostat=ios) xxx, elt
    if (ios /= 0) call error('read (str), #'//trim(string(ios)))
    call insert_sorted(el_type, elt, nelt, fit=.false.) !at the end, nelt is the total number of element types
  end do
  !allocate pmh%el(1)%el
  allocate(m%el(nelt), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_gmsh/load_gmsh) Unable to allocate piece: '//trim(cad))
  do i = 1, nelt
      m%el(i)%type = id4gmsh(el_type(i))
    if (id4gmsh(el_type(i)) == 0) call error('(module_gmsh/load_gmsh) unsupported GMSH element type: '//&
    trim(gmsh_elm_type(el_type(i))))
    if (FEDB(id4gmsh(i))%lnn >= size(tmp,1)) call error('(module_gmsh/load_gmsh) temporary array tmp cannot store nodes; '//&
    &'please, increase the tmp dimension and compile again.')
  end do
  !invert el_type: inv_el_type
  call alloc(inv_el_type, maxval(el_type))
  do i = 1, nelt
    inv_el_type(el_type(i)) = i
  end do
  !save element info
  rewind(iu)
  res = search_mark(iu, '$Elements')
  read (unit=iu, fmt=*, iostat=ios) nel
  tmp = 0
  do k = 1, nel
    read (unit=iu, fmt=*, iostat=ios) xxx, elt, ntags, (tag(i), i=1,ntags), tmp(1:FEDB(id4gmsh(elt))%lnn)
    if (ntags < 2) call error('(module_gmsh/load_gmsh) Gmsh requires at least two tags per element, failed at element #'//&
    &trim(string(k)))
    i = inv_el_type(elt) !index of m%el() where elt was saved
    tmp(1:FEDB(m%el(i)%type)%lnn) = new4old(tmp(1:FEDB(m%el(i)%type)%lnn)) !use the new node indices
    tmp(FEDB(m%el(i)%type)%lnn+1) = tag(1) !save the (first found) physical tag as reference, initially in nn
    if (find_sorted(2, m%el(i)%nn, tmp(1:FEDB(m%el(i)%type)%lnn), m%el(i)%nel) <= 0) then !this element was not previously saved
      call insert_sorted(2, m%el(i)%nn, tmp(1:FEDB(m%el(i)%type)%lnn+1), m%el(i)%nel, fit=[.true., .false.])
    end if
  end do
  !save references in el()%ref and reduce el()%nn
  do i = 1, size(m%el, 1)
    call set(m%el(i)%ref, m%el(i)%nn(FEDB(m%el(i)%type)%lnn+1,1:m%el(i)%nel), [(k, k=1,m%el(i)%nel)], fit=.true.)
    call reduce(m%el(i)%nn, FEDB(m%el(i)%type)%lnn, m%el(i)%nel)
  end do
end associate
close(iu)
call build_vertices(pmh)
end subroutine

!-----------------------------------------------------------------------
! save_gmsh: save a Gmsh MSH file from a PMH structure
!-----------------------------------------------------------------------
subroutine save_gmsh(outfile, iu, pmh)
character(*),   intent(in)    :: outfile
integer,        intent(in)    :: iu
type(pmh_mesh), intent(inout) :: pmh

integer :: i, ipp, ip, ig, k, j, id4pmh(16), valid_fe(12), res, ios, nel, nnod, prev_nel
integer, allocatable :: piece2save(:), nel_piece(:), nnod_piece(:)
real(real64), allocatable :: znod(:,:)
character(maxpath) :: str, cad
logical :: all_P1

!id4pmh("pmh element id") = "gmsh element id"
id4pmh(check_fe(.true.,   1, 1,  0, 0)) = 15 !Node
id4pmh(check_fe(.true.,   2, 2,  1, 0)) =  1 !Edge, Lagrange P1
id4pmh(check_fe(.false.,  3, 2,  1, 0)) =  8 !Edge, Lagrange P2
id4pmh(check_fe(.true.,   3, 3,  3, 0)) =  2 !Triangle, Lagrange P1
id4pmh(check_fe(.false.,  6, 3,  3, 0)) =  9 !Triangle, Lagrange P2
id4pmh(check_fe(.true.,   4, 4,  4, 0)) =  3 !Quadrangle, Lagrange P1
id4pmh(check_fe(.false.,  8, 4,  4, 0)) = 16 !Quadrangle, Lagrange P2
id4pmh(check_fe(.true.,   4, 4,  6, 4)) =  4 !Tetrahedron, Lagrange P1
id4pmh(check_fe(.false., 10, 4,  6, 4)) = 11 !Tetrahedron, Lagrange P2
id4pmh(check_fe(.true.,   8, 8, 12, 6)) =  5 !Hexahedron, Lagrange P1
id4pmh(check_fe(.false., 20, 8, 12, 6)) = 17 !Hexahedron, Lagrange P2
id4pmh(check_fe(.true.,   6, 6,  9, 5)) =  6 !Wedge, Lagrange P1

!valid elements types to save a MFM mesh (all but RT, ND)
valid_fe = [check_fe(.true.,   1, 1,  0, 0), & !Node
            check_fe(.true.,   2, 2,  1, 0), & !Edge, Lagrange P1
            check_fe(.false.,  3, 2,  1, 0), & !Edge, Lagrange P2
            check_fe(.true.,   3, 3,  3, 0), & !Triangle, Lagrange P1
            check_fe(.false.,  6, 3,  3, 0), & !Triangle, Lagrange P2
            check_fe(.true.,   4, 4,  4, 0), & !Quadrangle, Lagrange P1
            check_fe(.false.,  8, 4,  4, 0), & !Quadrangle, Lagrange P2
            check_fe(.true.,   4, 4,  6, 4), & !Tetrahedron, Lagrange P1
            check_fe(.false., 10, 4,  6, 4), & !Tetrahedron, Lagrange P2
            check_fe(.true.,   8, 8, 12, 6), & !Hexahedron, Lagrange P1
            check_fe(.false., 20, 8, 12, 6), & !Hexahedron, Lagrange P2
            check_fe(.true.,   6, 6,  9, 5)]   !Wedge, Lagrange P1

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

!testing
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_gmsh/save_gmsh) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (find_first(valid_fe, tp) == 0) then
        call info('(module_gmsh/save_gmsh) element type '//trim(FEDB(tp)%desc)//' found; those elements cannot be saved'//&
        &' in Gmsh format and they will be discarded')
        cycle
      end if
    end associate
  end do
end do

!calculate number of elements (nel) and number of nodes (nnod) for selected pieces
if (allocated(nel_piece)) deallocate(nel_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(nel_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_gmsh/save_gmsh) Unable to allocate variable nel_piece: '//trim(cad))
if (allocated(nnod_piece)) deallocate(nnod_piece); !nnod_piece(ipp): global numbering for the last node  of piece #ipp
allocate(nnod_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_gmsh/save_gmsh) Unable to allocate variable nnod_piece: '//trim(cad))
nel_piece(0) = 0; nnod_piece(0) = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  nnod_piece(ipp) = nnod_piece(ipp-1) + pmh%pc(ip)%nnod
  nel_piece(ipp)  =  nel_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    nel_piece(ipp) =  nel_piece(ipp) + pmh%pc(ip)%el(ig)%nel
  end do
end do
nnod = nnod_piece(size(piece2save,1))
nel  =  nel_piece(size(piece2save,1))

open (unit=iu, file=outfile, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('save/open, #'//trim(string(ios)))

!store header
write(iu, '(a)') '$MeshFormat'
write(iu, '(a)') '2.2 0 8'
write(iu, '(a)') '$EndMeshFormat'

!store nodes
write(iu, '(a)') '$Nodes'
write(iu, *)      nnod
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  call build_node_coordinates(pmh%pc(ip), ip, all_P1, znod)
  if (.not. all_P1) then
    do j = 1, pmh%pc(ip)%nnod
      write(iu, *) nnod_piece(ipp-1)+j, (znod(i,j), i = 1,pmh%pc(ip)%dim), (0._real64, i = 1,3-pmh%pc(ip)%dim)
    end do
  else
    do j = 1, pmh%pc(ip)%nver
      write(iu, *) nnod_piece(ipp-1)+j, (pmh%pc(ip)%z(i,j), i = 1,pmh%pc(ip)%dim), (0._real64, i = 1,3-pmh%pc(ip)%dim)
    end do
  end if
end do
write(iu, '(a)') '$EndNodes'
!store elements
write(iu, '(a)') '$Elements'
write(iu, *)      nel
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  prev_nel = nel_piece(ipp-1)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (.not. FEDB(elg%type)%nver_eq_nnod) then
        do k = 1, elg%nel
          write(iu, *) prev_nel+k, id4pmh(elg%type), 2, elg%ref(k), prev_nel+k, (nnod_piece(ipp-1)+elg%nn(i,k), &
          i = 1,FEDB(elg%type)%lnn)
        end do
      else
        do k = 1, elg%nel
          write(iu, *) prev_nel+k, id4pmh(elg%type), 2, elg%ref(k), prev_nel+k, (nnod_piece(ipp-1)+elg%mm(i,k), &
          i = 1,FEDB(elg%type)%lnv)
        end do
      end if
      prev_nel = prev_nel + elg%nel
    end associate
  end do
end do
write(iu, '(a)') '$EndElements'
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
