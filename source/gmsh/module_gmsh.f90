module module_gmsh
!-----------------------------------------------------------------------
! Module to manage Gmsh MSH meshes
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 26/05/2014
!
! PUBLIC PROCEDURES:
! load_gmsh: load a Gmsh ASCII MSH  file into a PMH structure
!-----------------------------------------------------------------------
use module_compiler_dependant, only: iostat_end
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_convers, only: string, lcase, adjustlt
use module_alloc, only: alloc, dealloc, insert_sorted, set, insert_col_sorted, find_col_sorted, reduce
use module_fe_database_pmh, only: FEDB, check_fe
use module_pmh, only: pmh_mesh, build_vertices
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
integer :: res, i, j, k, ios, nel, elt, nelt, id4gmsh(93), xxx, ntags, tag(100), tmp(100)
integer, allocatable :: el_type(:), inv_el_type(:), new4old(:), old4new(:)
character(maxpath) :: cad

!id4gmsh("gmsh element id") = "pmh element id"
id4gmsh(15) = check_fe(.true.,      1, 1,  0, 0) ! 1
id4gmsh( 1) = check_fe(.true.,      2, 2,  1, 0) ! 2 
id4gmsh( 8) = check_fe(.false.,     3, 2,  1, 0) ! 3
id4gmsh( 2) = check_fe(.true.,      3, 3,  3, 0) ! 4
id4gmsh( 9) = check_fe(.false.,     6, 3,  3, 0) ! 5
!id4gmsh( ) = check_fe(.false.,     3, 3,  3, 0) ! 6 (tri rt)
id4gmsh( 3) = check_fe(.true.,      4, 4,  4, 0) ! 7
id4gmsh(16) = check_fe(.false.,     8, 4,  4, 0) ! 8
id4gmsh( 4) = check_fe(.true.,      4, 4,  6, 4) ! 9
id4gmsh(11) = check_fe(.false.,    10, 4,  6, 4) !10
!id4gmsh( ) = check_fe(.false.,     4, 4,  6, 4) !11 (tet rt)
!id4gmsh( ) = check_fe(.false.,     6, 4,  6, 4) !12 (tet nd)
id4gmsh( 5) = check_fe(.true.,      8, 8, 12, 6) !13
id4gmsh(17) = check_fe(.false.,    20, 8, 12, 6) !14
id4gmsh( 6) = check_fe(.true.,      6, 6,  9, 5) !15

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
    if (FEDB(id4gmsh(elt))%lnn >= size(tmp,1)) call error('(module_gmsh/load_gmsh) temporary array tmp cannot store nodes; '//&
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
    if (find_col_sorted(m%el(i)%nn, tmp(1:FEDB(m%el(i)%type)%lnn), m%el(i)%nel) <= 0) then !this element was not previously saved
      call insert_col_sorted(m%el(i)%nn, tmp(1:FEDB(m%el(i)%type)%lnn+1), m%el(i)%nel, fit=[.true., .false.])
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
