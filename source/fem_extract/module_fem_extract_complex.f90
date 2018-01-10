module module_fem_extract_complex_fcnv
!-----------------------------------------------------------------------
! Module for fem extractions
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Rodrigo Valina; Francisco Pena, fran.pena@usc.es
! Last update: 17/03/2011
!
! PUBLIC TYPE-BOUND PROCEDURES:
!   extract_mesh: extract a submesh from a global mesh selecting 
!   only the elements with some particular subdomain numbers
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
implicit none

!Constants
integer, parameter, private :: real64 = selected_real_kind(15, 307)  ! double precision format

!Private procedures
private :: extract_subdomain, extract_z, extract_field_v, extract_field_1, extract_field_n, extract_field_c
private :: extract_field_complex, cell2node_complex

interface extract_field; module procedure extract_field_complex; end interface
interface cell2node;     module procedure cell2node_complex;     end interface

contains

!--------------------------------------------------------------------
! extract_field_complex: extract a subfield from a global field
!--------------------------------------------------------------------
subroutine extract_field_complex(globv, globel, v, subv, type_, ncomp)
integer,         dimension(:), intent(in) :: globv  ! return global vertex index
integer,         dimension(:), intent(in) :: globel ! return global elem.  index
complex(real64), dimension(:), intent(in) :: v      ! field values (global mesh)
complex(real64), dimension(:), allocatable, intent(out) :: subv  ! field values (submesh)
character(len=*),              intent(in) :: type_  ! field type (node or cell)
integer, optional,             intent(in) :: ncomp  ! number of field components (1 by default) 
integer :: lncomp

lncomp = 1
if (present(ncomp)) lncomp = ncomp

if (trim(adjustl(type_)) == 'node') then
  call extract_field_n(v, subv, globv, lncomp)
elseif (trim(adjustl(type_)) == 'cell') then
  call extract_field_n(v, subv, globel, lncomp)
else
  print*, '(extract_field) field type: '//trim(type_)
  stop 'Not implemented'
endif
end subroutine

!--------------------------------------------------------------------
! cell2node: calculate a node field form a cell one
!--------------------------------------------------------------------
subroutine cell2node_complex(nver, mm, vc, vn)
integer,                         intent(in) :: nver ! total number of vertices
integer,         dimension(:,:), intent(in) :: mm   ! conectivities 
complex(real64), dimension(:),   intent(in) :: vc   ! cell field
complex(real64), dimension(:), allocatable, intent(out) :: vn ! node field
integer, dimension(nver) :: cell4node
integer :: ncomp, lnn, nel, i, j, k, l

lnn   = size(mm, 1)
nel   = size(mm, 2)
ncomp = size(vc)/nel
cell4node = 0
if (allocated(vn)) deallocate(vn); allocate(vn(ncomp*nver))
vn = (0._real64,0._real64)
do k = 1, nel
  do j = 1, lnn
    i = mm(j,k)
    cell4node(i) = cell4node(i) + 1
    do l = 1, ncomp
      vn(ncomp*(i-1)+l) = vn(ncomp*(i-1)+l) + vc(ncomp*(k-1)+l)
    enddo
  enddo
enddo
do i = 1, nver
  do l = 1, ncomp
    vn(ncomp*(i-1)+l) = vn(ncomp*(i-1)+l)/cell4node(i)
  enddo
enddo
end subroutine

!--------------------------------------------------------------------
! PRIVATE PROCEDURES
!--------------------------------------------------------------------
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
 ! marca vertices usados
 do i = 1, size(mm,2)
   do k = 1, size(subnum,1) ! modified to allow several nsd0
     if (sub(i) == subnum(k)) then
        noc = noc + 1
        temp(mm(:, i)) = .true.
        exit
     endif
   enddo
 enddo

 nov = 0
 ! conta numero de vertice usados
 do i = 1, size(temp,1)
    if (temp(i)) nov = nov + 1
 enddo

 allocate(mm_new(size(mm,1),noc))
 allocate(map(nov))
 allocate(mapel(noc))
 j = 1
 ! calcula novo mm (falta cambiar indices)
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
 ! calcula punteiros de vertice novo a vello e de vello a novo
 do i = 1 , vnum
    if (temp(i)) then ! se se utiliza...
        map(j) = i
        map_inv(i) = j
        j = j + 1
    endif
 enddo

 ! reescribir novo mm, cabiando indices
 do i = 1, size(mm_new,2)
    do j = 1, size(mm_new,1)
        mm_new(j,i) = map_inv(mm_new(j,i))
    enddo
 enddo
end subroutine

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

!--------------------------------------------------------------------
! pass map for point fields and mapel for cell fields
subroutine extract_field_v(f, f_new, map)
 real(real64), dimension(:,:), intent(in) :: f
 real(real64), dimension(:,:), allocatable, intent(out) :: f_new
 integer, dimension(:), intent(in) :: map
 integer :: i
    allocate(f_new(size(f,1),size(map,1)))
    do i = 1 , size(map,1)
        f_new(:,i) = f(:,map(i))
    enddo
end subroutine

!--------------------------------------------------------------------
! pass map for point fields and mapel for cell fields
subroutine extract_field_1(f, f_new, map)
 real(real64), dimension(:), intent(in) :: f
 real(real64), dimension(:), allocatable, intent(out) :: f_new
 integer, dimension(:), intent(in) :: map
 integer :: i
    allocate(f_new(size(map,1)))
    do i = 1 , size(map,1)
        f_new(i) = f(map(i))
    enddo
end subroutine

!--------------------------------------------------------------------
! pass map for point fields and mapel for cell fields
subroutine extract_field_n(f, f_new, map, n)
 complex(real64), dimension(:), intent(in) :: f
 complex(real64), dimension(:), allocatable, intent(out) :: f_new
 integer, dimension(:), intent(in) :: map
 integer, intent(in) :: n
 integer :: i, j
    allocate(f_new(size(map,1)*n))
    do i = 1 , size(map,1)
        do j = 1 , n
            f_new((i-1)*n+j) = f((map(i)-1)*n+j)
        enddo
    enddo
end subroutine

!--------------------------------------------------------------------
! pass map for point fields and mapel for cell fields
subroutine extract_field_c(f, f_new, map)
 real(real64), dimension(:,:), intent(in) :: f
 real(real64), dimension(:), allocatable, intent(out) :: f_new
 integer, dimension(:), intent(in) :: map
 integer :: i, j
    allocate(f_new(size(map,1)*size(f,1)))
    do i = 1 , size(map,1)
        do j = 1 , size(f,1)
            f_new((i-1)*size(f,1)+j) = f(j, map(i))
        enddo
    enddo
end subroutine

end module
