module module_fem_extract_fcnv
!-----------------------------------------------------------------------
! Module for fem extractions
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Rodrigo Valina; Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 17/05/2014
!
! PUBLIC TYPE-BOUND PROCEDURES:
!   extract_mesh: extract a submesh from a global mesh selecting only the elements
!                 with some particular subdomain numbers
!   extract_ref: extract references
!   extract_field: extract a subfield from a global field
!   cell2node: calculate a node field form a cell one
!
! USAGE:
!   1) EXTRACT MESH: 
!      call extract_mesh(nver, mm, z, nsd, nsd0, submm, subz, globv, globel)
!
!      a) The global mesh is given by nver, mm, z, nsd
!      b) The subd. numbers to choose are given by nsd0
!      c) The submesh is returned in submm, subz
!      d) The index vectors to know the old indices of vertices and elements
!      are returned in globv, globel
!
!      Additionally, references can be extracted after extracting mesh:
!      call extract_ref(nrv, nra, nrc, nsd, subnrv, subnra, subnrc, subnsd, globel)
!
!      a) The global  references are given by nrv, nra, nrc, nsd
!      c) The submesh references are returned by subnrv, subnra, subnrc, subnsd
!      d) The index vector to know the old indices of elements is given in globel
!
!   2) EXTRACT FIELD:
!      call extract_field(globv, globel, v, subv, type, ncomp)
!
!      a) The index vector to know the old indices of vertices is given by globv
!      b) The index vector to know the old indices of elements is given by globel
!      c) The global field is given by v
!      d) The subfield is returned in subv
!      d) The character argument type must be 'node' or 'cell'
!      d) The optional integer argument ncomp is the component number; set to 1 by default.
!
!      Examples:
!      call extract_field(globv, globel, temp,  sub_temp,  'node')    ! scalar    node field
!      call extract_field(globv, globel, veloc, sub_veloc, 'node', 3) ! 3D vector node field
!      call extract_field(globv, globel, joule, sub_joule, 'cell')    ! scalar    cell field
!      call extract_field(globv, globel, gradV, sub_gradV, 'cell', 3) ! 3D vector cell field
!
!   2) CELL to NODE:
!      call cell2node(nver, mm, vc, vn)
!
!      a) The mesh is given by nver, mm
!      a) The cell field is given by vn
!      b) The node field is returned in vc
!-----------------------------------------------------------------------
use module_fem_extract_real_fcnv
use module_fem_extract_complex_fcnv
implicit none

!Constants
integer, parameter, private :: real64 = selected_real_kind(15, 307)

!Private procedures
private :: extract_z

contains

!--------------------------------------------------------------------
! extract_ref: extract references
!--------------------------------------------------------------------
subroutine extract_ref(nrv, nra, nrc, nsd, subnrv, subnra, subnrc, subnsd, globel)
integer, allocatable, intent(in)  ::    nrv(:,:),    nra(:,:),    nrc(:,:),    nsd(:)
integer, allocatable, intent(out) :: subnrv(:,:), subnra(:,:), subnrc(:,:), subnsd(:)
integer,              intent(in)  :: globel(:)
integer :: k

!nrv
if (allocated(nrv)) then
  allocate(subnrv(size(nrv,1), size(globel,1)))
  do k = 1, size(globel,1)
    subnrv(:,k) = nrv(:,globel(k))
  enddo
end if
!nra
if (allocated(nra)) then
  allocate(subnra(size(nra,1), size(globel,1)))
  do k = 1, size(globel,1)
    subnra(:,k) = nra(:,globel(k))
  enddo
end if
!nrc
if (allocated(nrc)) then
  allocate(subnrc(size(nrc,1), size(globel,1)))
  do k = 1, size(globel,1)
    subnrc(:,k) = nrc(:,globel(k))
  enddo
end if
!nsd
if (allocated(nsd)) then
  allocate(subnsd(size(globel,1)))
  do k = 1, size(globel,1)
    subnsd(k) = nsd(globel(k))
  enddo
end if
end subroutine

!--------------------------------------------------------------------
! PRIVATE PROCEDURES
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! extract_subdomain: extract the submesh connectivities
!--------------------------------------------------------------------
subroutine extract_subdomain(vnum, mm, sub, subnum, mm_new, map, mapel)
integer, intent(in) :: vnum
integer, dimension(:,:), intent(in) :: mm
integer, dimension(:), intent(in) :: sub
integer, intent(in) :: subnum(:)
integer, dimension(:,:), allocatable, intent(out) :: mm_new
integer, dimension(:), allocatable, intent(out) :: map ! new vertex -> old vertex
integer, dimension(:), allocatable, intent(out) :: mapel ! new element -> old element
integer, dimension(vnum) :: map_inv ! old vertex -> new vertex
logical, dimension(vnum) :: temp
integer :: i, j, nov, noc, k

temp = .false.
noc = 0
!marca vertices usados
do i = 1, size(mm,2)
  do k = 1, size(subnum,1) !modified to allow several nsd0
    if (sub(i) == subnum(k)) then
       noc = noc + 1
       temp(mm(:, i)) = .true.
       exit
    endif
  enddo
enddo

!conta numero de vertice usados
nov = 0
do i = 1, size(temp,1)
  if (temp(i)) nov = nov + 1
enddo
allocate(mm_new(size(mm,1),noc))
allocate(map(nov))
allocate(mapel(noc))
j = 1
!calcula novo mm (falta cambiar indices)
do i = 1, size(mm,2)
  do k = 1, size(subnum,1) ! modified to allow several nsd0
    if (sub(i) == subnum(k)) then
       mapel(j) = i
       mm_new(:, j) = mm(:, i)
       j = j + 1
       exit 
   endif
  enddo
enddo
j = 1
!calcula punteiros de vertice novo a vello e de vello a novo
do i = 1 , vnum
  if (temp(i)) then ! se se utiliza...
    map(j) = i
    map_inv(i) = j
    j = j + 1
  endif
enddo
!reescribir novo mm, cabiando indices
do i = 1, size(mm_new,2)
  do j = 1, size(mm_new,1)
    mm_new(j,i) = map_inv(mm_new(j,i))
  enddo
enddo
end subroutine

!--------------------------------------------------------------------
! extract_z: extract vertex coordinates
!--------------------------------------------------------------------
subroutine extract_z(z, z_new, map)
real(real64), dimension(:,:), intent(in) :: z
real(real64), dimension(:,:), allocatable, intent(out) :: z_new
integer, dimension(:), intent(in) :: map
integer :: i

allocate(z_new(size(z,1),size(map,1)))
do i = 1 , size(map,1)
  z_new(:,i) = z(:,map(i))
enddo
end subroutine

end module
