module module_dataset_2412
!-----------------------------------------------------------------------
! Module for dataset_2411 class
! Last update: 3/04/2010
!
! Name:   Elements
! Record 1:        FORMAT(6I10)
!                  Field 1       -- element label
!                  Field 2       -- fe descriptor id
!                  Field 3       -- physical property table number
!                  Field 4       -- material property table number
!                  Field 5       -- color
!                  Field 6       -- number of nodes on element
!
! Record 2:  *** FOR NON-BEAM ELEMENTS ***
!                  FORMAT(8I10)
!                  Fields 1-n    -- node labels defining element
!
! Record 2:  *** FOR BEAM ELEMENTS ONLY ***
!                  FORMAT(3I10)
!                  Field 1       -- beam orientation node number
!                  Field 2       -- beam fore-end cross section number
!                  Field 3       -- beam  aft-end cross section number
!
! Record 3:  *** FOR BEAM ELEMENTS ONLY ***
!                  FORMAT(8I10)
!                  Fields 1-n    -- node labels defining element
!
! Records 1 and 2 are repeated for each non-beam element in the model.
! Records 1 - 3 are repeated for each beam element in the model.
!
! Example:
!
!     -1
!   2412
!          1        11         1      5380         7         2
!          0         1         1
!          1         2
!          2        21         2      5380         7         2
!          0         1         1
!          3         4
!          3        22         3      5380         7         2
!          0         1         2
!          5         6
!          6        91         6      5380         7         3
!         11        18        12
!          9        95         6      5380         7         8
!         22        25        29        30        31        26        24        23
!         14       136         8         0         7         2
!         53        54
!         36       116        16      5380         7        20
!        152       159       168       167       166       158       150       151
!        154       170       169       153       157       161       173       172
!        171       160       155       156
!     -1
!-----------------------------------------------------------------------
use module_ALLOC
use module_dataset
use module_mesh_unv
use module_FE_DB
use module_cells
implicit none

!Private procedures
private :: check_and_set, reorder_nodes, find_prv, swap

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! read: read dataset 2412
!-----------------------------------------------------------------------
subroutine read_2412(iu, m, maxdim, is_opt)

  integer,                                intent(in)    :: iu     !unit number for unvfile
  type(mfm_mesh), allocatable, dimension(:), intent(inout) :: m      !meshes
  integer,                                intent(inout) :: maxdim !maximum founded dimension
  logical,                                intent(in)    :: is_opt !-is option
  integer :: ios, Field1, Field2, F3, F4, F5, Field6, F1, F2, i, d
  integer, allocatable, dimension(:) :: nn !node numbers
  logical :: fit(2)

  maxdim = 0
  call FE_DB_init() !initialize FE database
  if(.not. allocated(m))  call error('dataset_2412/read, mesh(es) not allocated')
  do
    if (is_dataset_delimiter(iu, back=.true.)) exit

!   Record 1
    read (unit=iu, fmt='(6I10)', iostat = ios) &
    Field1, & !element label
    Field2, & !fe descriptor id
    F3,     & !physical property table number
    F4,     & !material property table number
    F5,     & !color
    Field6    !number of nodes on element
    if (ios /= 0) call error('dataset_2412/read, #'//trim(string(ios)))

!   check whether the FE is in the database
    if (len_trim(FE_DB(Field2)%desc)==0) &
    call error('dataset_2412/read, fe descriptor id '//trim(string(Field2))//' not implemented')
    d = FE_DB(Field2)%DIM
    maxdim = max(maxdim, d)
!   set local quantities (an error occurs if they differ from values already stored)
    call check_and_set('LNN', m(d)%LNN, FE_DB(Field2)%LNN)
    call check_and_set('LNV', m(d)%LNV, FE_DB(Field2)%LNV)
    call check_and_set('LNE', m(d)%LNE, FE_DB(Field2)%LNE)
    call check_and_set('LNF', m(d)%LNF, FE_DB(Field2)%LNF)
    call set_subelements(m(d)) !determine sub-elements of the mesh
    m(d)%FEtype = FE_DB(Field2)%desc

!   Record 2:  *** FOR BEAM ELEMENTS ONLY ***
    if (FE_DB(Field2)%Beam) then
      read (unit=iu, fmt='(3I10)', iostat = ios) &
      F1, & !beam orientation node number
      F2, & !beam fore-end cross section number
      F3    !beam  aft-end cross section number
      if (ios /= 0) call error('dataset_2412/read, #'//trim(string(ios)))
    end if

!   Record 2 (FOR NON-BEAM ELEMENTS) or 3 (FOR BEAM ELEMENTS ONLY)
    call alloc(nn, m(d)%LNN)
    read (unit=iu, fmt='(8I10)', iostat = ios) &
    (nn(i), i = 1, m(d)%LNN) !Fields 1-n    -- node labels defining element
    if (ios /= 0) call error('dataset_2412/read, #'//trim(string(ios)))

!   reorder nodes to obtain vertices first and counter-clockwise orientation
    call reorder_nodes(m(size(m,1)), m(d), nn, is_opt)

!   copy nodes to mesh (add an element after the last one)
    fit = [.true., .false.]
    call set(2, m(d)%id, nn, m(d)%nl+1, fit)
    m(d)%nl = m(d)%nl + 1

!   set the new element in fes list of elements
    fit = [.true., .false.]
    call set(2, cells, [d, m(d)%nl], Field1, fit)
  end do
  call reduce(m(maxdim)%id, m(maxdim)%LNN, m(maxdim)%nl)
end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! check_and_set: check whether x and y differs; if not assigns y to x
!-----------------------------------------------------------------------
subroutine check_and_set(str, x, y)
  character(len=*), intent(in) :: str
  integer,          intent(inout) :: x
  integer,          intent(in)    :: y

!  if (y == 0) then
!    call error('dataset_2412/read/check_and_set, fe has null '//str//' in database')
  if (x>0 .and. x/=y) then
    call error('dataset_2412/read/check_and_set, '//str//' mesh differs from Field6')
  else
    x = y
  end if

end subroutine

!-----------------------------------------------------------------------
! reorder_nodes: reorder nn to obtain vertices first and counter-clockwise orientation
! REMARK: only valid for non-isoparametric triangles (P1 & P2) and tetrahedra (P1 & P2)
!-----------------------------------------------------------------------
subroutine reorder_nodes(mx, m, nn, is_opt)
type(mfm_mesh), intent(in)    :: mx     !mesh with coordinates
type(mfm_mesh), intent(in)    :: m      !mesh
integer,     intent(inout) :: nn(:)  !nodes
logical,     intent(in)    :: is_opt !-is option
integer :: i, j, l, vf, newnn(10)
real(real64) :: a2(3), a3(3), a4(3)
! Init newnn
newnn=-1

if (is_opt) then
  !Input file contains P2 isoparametrical elements that must be read using Salome ordering (vertices and midpoints sandwiched)
  i = 1
  do j = 1, m%LNN, 2; newnn(i) = nn(j); i = i+1; end do
  do j = 2, m%LNN, 2; newnn(i) = nn(j); i = i+1; end do
  nn = newnn(1:size(nn,1))
elseif (m%LNN /= m%LNV) then
  !identify vertices (non midpoints)
  vf = 0
  NODES: do i = 1, m%LNN !nodes
    do j = 1, m%LNN !first possible vertex
      if (j /= i) then
        do l = j+1, m%LNN !second possible vertex
          if (l /= i) then
            if ( maxval(abs((mx%xd(:,nn(j))+mx%xd(:,nn(l)))/2-mx%xd(:,nn(i)))) < 1e3*epsilon(mx%xd) ) cycle NODES
            !the check: maxval(abs(mx%xd(:,nn(j))-mx%xd(:,nn(l)))) > 1e3*epsilon(mx%xd), to see whether (j,l) is a singular egde
            ! is not made anymore; instead of that, when the number of non mid-points is greater then LNV, an error arises to warn
            ! that P2 elements are isoparametrical
          end if
        end do
      end if
    end do
    !if it is not a midpoint, it is a vertex
    vf = vf + 1
    newnn(vf) = nn(i)
    !if (vf >= m%LNV) exit NODES !found enough vertices (this is useful to deal with singular edges)
    if (vf > m%LNV)  call error('dataset_2412/reorder_nodes, too many non midpoints in an element (isoparametric P2 elements &
    &are not supported)')
  end do NODES
  !identify midpoints
  do i = 1, m%LNN !nodes
    if (find_first(newnn == nn(i))>0) cycle !it is a vertex
    do j = 1, size(m%edge,2) !avoid the use of LNE: it is 0 for beams
      if ( maxval(abs((mx%xd(:,newnn(m%edge(1,j)))+mx%xd(:,newnn(m%edge(2,j))))/2-mx%xd(:,nn(i)))) < 1e3*epsilon(mx%xd) ) then
        newnn(m%LNV+j) = nn(i)
        exit
      end if
    end do
  end do
  nn = newnn(1:size(nn,1))
end if

!ensure counter-clockwise orientation
if (m%DIM == 2 .and. m%LNV == 3 .and. m%LNE == 3) then !triangles
  a2 = mx%xd(:,nn(2))-mx%xd(:,nn(1))
  a3 = mx%xd(:,nn(3))-mx%xd(:,nn(1))
!  if (abs(mx%xd(3,nn(1))) > 1e3*epsilon(mx%xd) .or. &
!      abs(mx%xd(3,nn(2))) > 1e3*epsilon(mx%xd) .or. &
!      abs(mx%xd(3,nn(3))) > 1e3*epsilon(mx%xd) .and. .not. tria_non_coplanar) tria_non_coplanar = .true.
  if ( a2(1)*a3(2)-a2(2)*a3(1) < 0 ) then
    call swap(nn(2), nn(3))
    if (m%LNN == 6) call swap(nn(4), nn(6)) !triangles P2
  end if
elseif (m%DIM == 3 .and. m%LNV == 4 .and. m%LNE == 6 .and. m%LNF == 4) then !tetrahedra
  a2 = mx%xd(:,nn(2))-mx%xd(:,nn(1))
  a3 = mx%xd(:,nn(3))-mx%xd(:,nn(1))
  a4 = mx%xd(:,nn(4))-mx%xd(:,nn(1))
  if ( a2(1)*a3(2)*a4(3) + a2(3)*a3(1)*a4(2) + a2(2)*a3(3)*a4(1) &
      -a2(3)*a3(2)*a4(1) - a2(2)*a3(1)*a4(3) - a2(1)*a3(3)*a4(2) < 0 ) then
    call swap(nn(2), nn(3))
    if (m%LNN == 10) then !tetrahedra P2
      call swap(nn(5),  nn(7))
      call swap(nn(9), nn(10))
    end if
  end if
end if

end subroutine

!-----------------------------------------------------------------------
! find (not  used): find the first occurrence of a value in a vector
!-----------------------------------------------------------------------
function find_prv(v) result(res)
logical, dimension(:), intent(in) :: v
integer :: res, tmp(size(v,1)), a(1)

tmp = 0
where (v); tmp = 1; end where
if (all(tmp==0)) then
  res = 0
else
  a = maxloc(tmp); res = a(1)
end if
end function

!-----------------------------------------------------------------------
! swap: swap the value of two variables
!-----------------------------------------------------------------------
subroutine swap(a,b)
integer :: a,b,c
c = a; a = b; b = c
end subroutine

end module
