module module_fem_extract_real_fcnv
!-----------------------------------------------------------------------
! Module for fem extractions (real64 fields)
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Rodrigo Valina; Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 17/05/2014
!
! PUBLIC TYPE-BOUND PROCEDURES:
!   extract_field: extract a subfield from a global field
!   cell2node: calculate a node field form a cell one
!-----------------------------------------------------------------------
implicit none

!Constants
integer, parameter, private :: real64 = selected_real_kind(15, 307)

!Private procedures
private :: extract_field_n, extract_field_real, cell2node_real
private :: extract_field_v, extract_field_1, extract_field_c !(unused procedures)

!Interfaces
interface extract_field; module procedure extract_field_real; end interface
interface cell2node;     module procedure cell2node_real;     end interface

contains

!--------------------------------------------------------------------
! extract_field: extract a subfield from a global field
!--------------------------------------------------------------------
subroutine extract_field_real(globv, globel, v, subv, type_, ncomp)
integer,      dimension(:), intent(in) :: globv  ! return global vertex index
integer,      dimension(:), intent(in) :: globel ! return global elem.  index
real(real64), dimension(:), intent(in) :: v      ! field values (global mesh)
real(real64), dimension(:), allocatable, intent(out) :: subv  ! field values (submesh)
character(len=*),           intent(in) :: type_  ! field type (node or cell)
integer, optional,          intent(in) :: ncomp  ! number of field components (1 by default) 
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
subroutine cell2node_real(nver, mm, vc, vn)
integer,                      intent(in) :: nver ! total number of vertices
integer,      dimension(:,:), intent(in) :: mm   ! conectivities 
real(real64), dimension(:),   intent(in) :: vc   ! cell field
real(real64), dimension(:), allocatable, intent(out) :: vn ! node field
integer, dimension(nver) :: cell4node
integer :: ncomp, lnn, nel, i, j, k, l

lnn   = size(mm, 1)
nel   = size(mm, 2)
ncomp = size(vc)/nel
cell4node = 0
if (allocated(vn)) deallocate(vn); allocate(vn(ncomp*nver))
vn = 0._real64
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
! extract_field_n: pass map for point fields and mapel for cell fields
!--------------------------------------------------------------------
subroutine extract_field_n(f, f_new, map, n)
 real(real64), dimension(:), intent(in) :: f
 real(real64), dimension(:), allocatable, intent(out) :: f_new
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
! UNUSED PROCEDURES
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! pass map for point fields and mapel for cell fields
!--------------------------------------------------------------------
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
!--------------------------------------------------------------------
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
!--------------------------------------------------------------------
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
