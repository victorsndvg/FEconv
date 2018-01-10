module module_dataset_2412_fcnv
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
use basicmod
use module_alloc_common_bmod, only: search_multiple !basicmod does not directly provides this procedure
use module_dataset_fcnv
use module_mesh_unv_fcnv
use module_FE_DB_fcnv
use module_cells_fcnv
use module_pmh_fcnv, only:piece
use module_fe_database_pmh_fcnv, only:FEDB, check_fe
implicit none

!Private procedures
private :: check_and_set, reorder_nodes, reorder_nodes_pmh, find_prv, swap

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! read: read dataset 2412
!-----------------------------------------------------------------------
subroutine read_2412(iu, pc, maxdim, els_loc)
integer,              intent(in)    :: iu           !unit number for unvfile
type(piece),          intent(inout) :: pc           !piece
integer,              intent(inout) :: maxdim       !maximum founded dimension
integer, allocatable, intent(out)   :: els_loc(:,:) !([elgroup, pos], element label)
integer :: ios, Field1, Field2, F3, F4, F5, Field6, F1, F2, i, d, gcounter
integer, allocatable :: nn(:) , elsingroup(:)
type(piece) :: auxpiece
logical :: fit(2)

call info('Reading mesh connectivities...')
maxdim = 0
gcounter = 1
if(.not. allocated(elsingroup)) allocate(elsingroup(size(FE_DB,1)))
elsingroup = 0
call FE_DB_init() !initialize FE database

do
  if (is_dataset_delimiter(iu, back=.true.)) exit
  !record 1
  read (unit=iu, fmt='(6I10)', iostat = ios) &
  Field1, & !element label
  Field2, & !fe descriptor id
  F3,     & !physical property table number
  F4,     & !material property table number
  F5,     & !color
  Field6    !number of nodes on element
  if (ios /= 0) call error('dataset_2412/read, #'//trim(string(ios)))
  !check whether the FE is in the database
  if (len_trim(FE_DB(Field2)%desc)==0) &
  call error('dataset_2412/read, fe descriptor id '//trim(string(Field2))//' not implemented')
  d = FE_DB(Field2)%DIM
  maxdim = max(maxdim, d)
  !assign element group by element type
  if(elsingroup(Field2) == 0) then
    elsingroup(Field2) = gcounter
    gcounter = gcounter + 1
  endif
  !allocate pc%el
  if(.not. allocated(pc%el)) then
    allocate(pc%el(elsingroup(Field2)))
    pc%el(elsingroup(Field2))%type = &
    & check_fe(FE_DB(Field2)%LNN==FE_DB(Field2)%LNV, FE_DB(Field2)%LNN, &
    & FE_DB(Field2)%LNV, FE_DB(Field2)%LNE, FE_DB(Field2)%LNF)
    call info('  Element type: '//trim(FEDB(pc%el(elsingroup(Field2))%type)%desc))
  endif
  !extend pc%el
  if(size(pc%el,1) < elsingroup(Field2)) then
    if(allocated(auxpiece%el)) deallocate(auxpiece%el)
    allocate(auxpiece%el(elsingroup(Field2)))
    auxpiece%el(1:size(pc%el,1)) = pc%el(:)
    call move_alloc(from=auxpiece%el,to=pc%el)
    pc%el(elsingroup(Field2))%type = &
    & check_fe(FE_DB(Field2)%LNN==FE_DB(Field2)%LNV, FE_DB(Field2)%LNN, &
    & FE_DB(Field2)%LNV, FE_DB(Field2)%LNE, FE_DB(Field2)%LNF)
    call info('  Element type: '//trim(FEDB(pc%el(elsingroup(Field2))%type)%desc))
  endif

  !record 2:  *** FOR BEAM ELEMENTS ONLY ***
  if (FE_DB(Field2)%Beam) then
    read (unit=iu, fmt='(3I10)', iostat = ios) &
    F1, & !beam orientation node number
    F2, & !beam fore-end cross section number
    F3    !beam  aft-end cross section number
    if (ios /= 0) call error('dataset_2412/read, #'//trim(string(ios)))
  end if
  !record 2 (FOR NON-BEAM ELEMENTS) or 3 (FOR BEAM ELEMENTS ONLY)
  call alloc(nn, FE_DB(Field2)%LNN)
  read (unit=iu, fmt='(8I10)', iostat = ios) &
  (nn(i), i = 1, FE_DB(Field2)%LNN) !Fields 1-n    -- node labels defining element
  !call reorder_nodes_pmh(pc%z, pc%el(elsingroup(Field2))%type, nn, is_opt)
  fit = [.true., .false.]
  pc%el(elsingroup(Field2))%nel = pc%el(elsingroup(Field2))%nel + 1
  call set(2, pc%el(elsingroup(Field2))%nn, nn, pc%el(elsingroup(Field2))%nel, fit)
  call set(2, els_loc, [elsingroup(Field2), pc%el(elsingroup(Field2))%nel], Field1, fit)
  if (ios /= 0) call error('dataset_2412/read, #'//trim(string(ios)))
end do

do i=1, size(pc%el,1)
  call reduce(pc%el(i)%nn, FEDB(pc%el(i)%type)%lnn, pc%el(i)%nel)
  allocate(pc%el(i)%ref(pc%el(i)%nel))
  pc%el(i)%ref = 0
enddo
end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! check_and_set (unused): check whether x and y differs; if not assigns y to x
!-----------------------------------------------------------------------
subroutine check_and_set(str, x, y)
character(len=*), intent(in) :: str
integer,          intent(inout) :: x
integer,          intent(in)    :: y

if (x>0 .and. x/=y) then
  call error('dataset_2412/read/check_and_set, '//str//' mesh differs from Field6')
else
  x = y
end if
end subroutine

!-----------------------------------------------------------------------
! reorder_nodes (unused): reorder nn to obtain vertices first and counter-clockwise orientation
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
            !is not made anymore; instead of that, when the number of non mid-points is greater then LNV, an error arises to warn
            !that P2 elements are isoparametrical
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
    if (find_first(newnn, nn(i))>0) cycle !it is a vertex
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
! reorder_nodes_pmh (unused): reorder nn to obtain vertices first and counter-clockwise orientation
! REMARK: only valid for non-isoparametric triangles (P1 & P2) and tetrahedra (P1 & P2)
!-----------------------------------------------------------------------
subroutine reorder_nodes_pmh(z, eltype, nn, is_opt)
real(real64),allocatable,dimension(:,:), intent(in) :: z     !mesh coordinates
integer,                            intent(in) :: eltype      !mesh
integer,                         intent(inout) :: nn(:)  !nodes
logical,                            intent(in) :: is_opt !-is option
integer, dimension(2,12)                       :: edge
integer :: i, j, l, vf, newnn(20), lnn,lnv,lne,lnf,dim
real(real64) :: a2(3), a3(3), a4(3)
! Init newnn
newnn=-1

lnn = FEDB(eltype)%lnn
lnv = FEDB(eltype)%lnv
lne = FEDB(eltype)%lne
lnf = FEDB(eltype)%lnf
dim = FEDB(eltype)%tdim
edge = FEDB(eltype)%edge

if (is_opt) then
  !Input file contains P2 isoparametrical elements that must be read using Salome ordering (vertices and midpoints sandwiched)
  i = 1
  do j = 1, lnn, 2; newnn(i) = nn(j); i = i+1; end do
  do j = 2, lnn, 2; newnn(i) = nn(j); i = i+1; end do
  nn = newnn(1:size(nn,1))
elseif (lnn /= lnv) then
  !identify vertices (non midpoints)
  vf = 0
  NODES: do i = 1, lnn !nodes
    do j = 1, lnn !first possible vertex
      if (j /= i) then
        do l = j+1, lnn !second possible vertex
          if (l /= i) then
            if ( maxval(abs((z(:,nn(j))+z(:,nn(l)))/2-z(:,nn(i)))) < 1e3*epsilon(0._real64) ) cycle NODES
            !the check: maxval(abs(z(:,nn(j))-z(:,nn(l)))) > 1e3*epsilon(0._real64), to see whether (j,l) is a singular egde
            !is not made anymore; instead of that, when the number of non mid-points is greater then LNV, an error arises to warn
            !that P2 elements are isoparametrical
          end if
        end do
      end if
    end do
    !if it is not a midpoint, it is a vertex
    vf = vf + 1
    newnn(vf) = nn(i)
    !if (vf >= lnv) exit NODES !found enough vertices (this is useful to deal with singular edges)
    if (vf > lnv)  call error('dataset_2412/reorder_nodes, too many non midpoints in an element (isoparametric P2 elements &
    &are not supported)')
  end do NODES
  !identify midpoints
  do i = 1, lnn !nodes
    if (find_first(newnn, nn(i))>0) cycle !it is a vertex
    do j = 1, lne
      if ( maxval(abs((z(:,newnn(edge(1,j)))+z(:,newnn(edge(2,j))))/2-z(:,nn(i)))) < 1e3*epsilon(0._real64) ) then
        newnn(lnv+j) = nn(i)
        exit
      end if
    end do
  end do
  nn = newnn(1:size(nn,1))
end if

!ensure counter-clockwise orientation
if (dim == 2 .and. lnv == 3 .and. lne == 3) then !triangles
  a2 = z(:,nn(2))-z(:,nn(1))
  a3 = z(:,nn(3))-z(:,nn(1))
!  if (abs(z(3,nn(1))) > 1e3*epsilon(0._real64) .or. &
!      abs(z(3,nn(2))) > 1e3*epsilon(0._real64) .or. &
!      abs(z(3,nn(3))) > 1e3*epsilon(0._real64) .and. .not. tria_non_coplanar) tria_non_coplanar = .true.
  if ( a2(1)*a3(2)-a2(2)*a3(1) < 0 ) then
    call swap(nn(2), nn(3))
    if (lnn == 6) call swap(nn(4), nn(6)) !triangles P2
  end if
elseif (dim == 3 .and. lnv == 4 .and. lne == 6 .and. lnf == 4) then !tetrahedra
  a2 = z(:,nn(2))-z(:,nn(1))
  a3 = z(:,nn(3))-z(:,nn(1))
  a4 = z(:,nn(4))-z(:,nn(1))
  if ( a2(1)*a3(2)*a4(3) + a2(3)*a3(1)*a4(2) + a2(2)*a3(3)*a4(1) &
      -a2(3)*a3(2)*a4(1) - a2(2)*a3(1)*a4(3) - a2(1)*a3(3)*a4(2) < 0 ) then
    call swap(nn(2), nn(3))
    if (lnn == 10) then !tetrahedra P2
      call swap(nn(5),  nn(7))
      call swap(nn(9), nn(10))
    end if
  end if
end if

end subroutine

!-----------------------------------------------------------------------
! find (unused): find the first occurrence of a value in a vector
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
