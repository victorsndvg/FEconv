module module_fem_extract_complex
!-----------------------------------------------------------------------
! Module for fem extractions (complex64 fields)
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
private :: extract_field_complex, cell2node_complex, extract_field_n

!Interfaces
interface extract_field; module procedure extract_field_complex; end interface
interface cell2node;     module procedure cell2node_complex;     end interface

contains

!--------------------------------------------------------------------
! extract_field: extract a subfield from a global field
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
! extract_field_n: pass map for point fields and mapel for cell fields
!--------------------------------------------------------------------
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

end module
