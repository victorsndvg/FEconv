module module_mesh_unv_fcnv
!-----------------------------------------------------------------------
! Module for mesh management
! Last update: 04/04/2010
!
! NOTES: 1) the mesh must be composed of a single finite element
!        2) the degrees of freedom are the vertices
!-----------------------------------------------------------------------
use basicmod, only: feed, empty, real64, maxpath, string, get_unit, error, alloc
implicit none

!Constants
!character(len = 255), private :: line = ' ' !line (formatted writing)

!Types
type mfm_mesh
! variables filled by data outside the Modulef mesh file
  integer :: DIM = 0 ! space dimension
  integer :: LNN = 0 ! local number of nodes
  integer :: LNV = 0 ! local number of vertices
  integer :: LNE = 0 ! local number of edges
  integer :: LNF = 0 ! local number of faces
  integer, dimension(:, :), allocatable :: edge  !local indices of vertices of the edges of an element (v. per edge x LNE)
  integer, dimension(:, :), allocatable :: face  !local indices of vertices of the faces of an element (v. per face x LNF)
  character(len=MAXPATH) :: filename = ' ' !file name
  integer                :: UNIT     = -1  !associated unit number
! variables filled by the Modulef mesh file
  character(len=MAXPATH)                     :: FEtype = ' ' !FE type
  integer                                    :: nl = 0 !global number of elements
  integer                                    :: nd = 0 !global number of elements
  integer                                    :: nv = 0 !global number of vertices
  integer,      dimension(:, :), allocatable :: id !nodes index array
  integer,      dimension(:, :), allocatable :: iv !vertices index array
  integer,      dimension(:, :), allocatable :: rv !vertices reference array
  integer,      dimension(:, :), allocatable :: re !edge reference array
  integer,      dimension(:, :), allocatable :: rf !face reference array
  real(real64), dimension(:, :), allocatable :: xv !vertices coordinates array
  real(real64), dimension(:, :), allocatable :: xd !nodes coordinates array
  integer,      dimension(:),    allocatable :: rl !element index array
end type


!Private procedures
private :: stringd !feed, empty, stringd

contains

!***********************************************************************
! MANAGEMENT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! dealloc: deallocate arrays
!-----------------------------------------------------------------------
subroutine dealloc_mesh(this)
type(mfm_mesh) :: this

if (allocated(this%id)) deallocate(this%id)
if (allocated(this%iv)) deallocate(this%iv)
if (allocated(this%rf)) deallocate(this%rf)
if (allocated(this%re)) deallocate(this%re)
if (allocated(this%rv)) deallocate(this%rv)
if (allocated(this%xv)) deallocate(this%xv)
if (allocated(this%xd)) deallocate(this%xd)
if (allocated(this%rl)) deallocate(this%rl)
this%nl = 0
this%nv = 0

end subroutine

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! set_subelements: sets which vertices define the edges and faces
!-----------------------------------------------------------------------
subroutine set_subelements(this)
type(mfm_mesh) :: this
logical :: stablished

stablished = .false.
if (this%DIM == 1) then !segments
  call alloc(this%edge, 2, 1)
  this%edge = reshape([1,2],[2,1])
  stablished = .true.
elseif (this%LNV == 3 .and. this%LNE == 3) then !triangles
  call alloc(this%edge, 2, 3)
  this%edge = reshape([1,2, 2,3, 3,1],[2,3])
  stablished = .true.
elseif (this%LNV == 4 .and. this%LNE == 4) then !quadrilateral
  call alloc(this%edge, 2, 4)
  this%edge = reshape([1,2, 2,3, 3,4, 4,1],[2,4])
  stablished = .true.
elseif (this%LNV == 4 .and. this%LNE == 6 .and. this%LNF == 4) then !tetrahedra
  call alloc(this%edge, 2, 6)
  call alloc(this%face, 3, 4)
  this%edge = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4],[2,6])
  this%face = reshape([1,3,2, 1,4,3, 1,2,4, 2,3,4],[3,4])
  stablished = .true.
else
  call error('mesh/set_subelements, this type of FE is not implemented')
end if
if (.not. stablished) call error('mesh/set_sub_elements, element type not recognized')

end subroutine

!-----------------------------------------------------------------------
! open: open mesh file
!-----------------------------------------------------------------------
subroutine open_mesh(this, filename)

type(mfm_mesh)                     :: this     !unv object
character(len=*), intent(in)    :: filename !unv file
integer :: ios

! open file
this%filename = filename
this%unit = get_unit()
open (unit=this%unit, file=this%filename, form='formatted', iostat=ios, position='rewind')
if (ios /= 0) call error('mesh/open, #'//trim(string(ios)))

end subroutine

!***********************************************************************
! OUTPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! save: save mesh
!-----------------------------------------------------------------------
subroutine write_mesh(this)
type(mfm_mesh) :: this
!integer, dimension(:,:), allocatable :: in !nodes index array, DISABLED
integer :: i, j, k, ln2, le2, lf2

!write data (nel, nnod, nver, dim, ...)
call feed(this%UNIT, string(this%nl));  call feed(this%UNIT, string(this%nd));   call feed(this%UNIT, string(this%nv))
call feed(this%UNIT, string(this%DIM)); call feed(this%UNIT, string(this%LNN));  call feed(this%UNIT, string(this%LNV))
call feed(this%UNIT, string(this%LNE)); call feed(this%UNIT, string(this%LNF));  call empty(this%UNIT)
!ln2: write nodes, only if (nnod /= nver)
ln2 = this%LNN; if (this%nd == this%nv) ln2 = 0
le2 = this%LNE; if (this%DIM < 2)       le2 = 0
lf2 = this%LNF; if (this%DIM < 3)       lf2 = 0
!arrays from file ([nn,if nnod/=nver], mm, [nrc,if dim==3], [nra,if dim==2], nrv, x)
do k = 1, this%nl;  do i = 1, ln2;      call feed(this%UNIT, string( this%id(i,k)));   end do; end do
do k = 1, this%nl;  do i = 1, this%LNV; call feed(this%UNIT, string( this%iv(i,k)));   end do; end do
do k = 1, this%nl;  do i = 1, lf2;      call feed(this%UNIT, string( this%rf(i,k)));   end do; end do
do k = 1, this%nl;  do i = 1, le2;      call feed(this%UNIT, string( this%re(i,k)));   end do; end do
do k = 1, this%nl;  do i = 1, this%LNV; call feed(this%UNIT, string( this%rv(i,k)));   end do; end do
do j = 1, this%nv;  do i = 1, this%DIM; call feed(this%UNIT, stringd(this%xv(i,j)));   end do; end do
call empty(this%UNIT)
!subdomain index array (nsd)
do k = 1, this%nl; call feed(this%UNIT, string(this%rl(k)));  end do
call empty(this%UNIT)
close(this%UNIT)

end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!--------------------------------------------------------------------
! save 'arg' in 'line'; if 'arg' is too large, write 'line' to file
!--------------------------------------------------------------------
!subroutine feed(iu, arg)
!
!  integer, intent(in) :: iu
!  character(len = *), intent(in) :: arg
!
!  if (len_trim(arg) > len(line) - len_trim(line) - 1) then
!    write(iu, '(a)') trim(line)
!    line = trim(adjustl(arg))
!  elseif (len_trim(line) == 0) then
!    line = trim(adjustl(arg))
!  else
!    line = trim(line)//' '//trim(adjustl(arg))
!  end if
!
!end subroutine

!--------------------------------------------------------------------
! empty 'line', writing it to file
!--------------------------------------------------------------------
!subroutine empty(iu)
!
!  integer, intent(in) :: iu
!
!  write(iu, '(a)') trim(line)
!  line = ' '
!
!end subroutine

!--------------------------------------------------------------------
! stringd: converts a numeric variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function stringd(i) result(res)

  real(real64), intent(in) :: i
  character(len=MAXPATH) :: res
  integer :: ios

  write(unit=res, fmt=*, iostat=ios) i
  if (ios /= 0) call error('string_dp/write')
  res = adjustl(res)

end function

end module
