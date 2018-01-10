module module_transform_fcnv
!-----------------------------------------------------------------------
! Module to transform Lagrange FE meshes to Lagrange P2, Raviart-Thomas
! and Whitney FE meshes
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 22/05/2013
!
! PUBLIC PROCEDURES:
!   lagr2l2: transforms a Lagrange FE mesh into a Lagrange P2 FE mesh
!   lagr2rt: transforms a Lagrange FE mesh into a Lagrange Raviart-Thomas (face) FE mesh
!   lagr2nd: transforms a Lagrange FE mesh into a Whitney (edge) FE mesh
!   to_l1:   transforms a mesh into a Lagrange P1 FE mesh
!-----------------------------------------------------------------------
use basicmod, only: maxpath, error, info, dealloc, alloc, sort, dealloc, alloc, insert_sorted, reduce, find_sorted
use module_vtu_fcnv, only: type_cell, edge_tetra, edge_tria, face_tetra
use module_cuthill_mckee_fcnv, only: bandwidth
use module_pmh_fcnv
implicit none

integer :: nparts = 0
integer, allocatable :: part(:,:) !part of a FE where new DOF is linked (face, edge,...)

contains

!-----------------------------------------------------------------------
! lagr2l2: transforms a Lagrange P1 FE mesh into a Lagrange P2 FE mesh
!-----------------------------------------------------------------------
subroutine lagr2l2(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
integer, intent(in)    :: nel        !global number of elements
integer, intent(inout) :: nnod       !global number of nodes
integer, intent(inout) :: nver       !global number of vertices
integer, intent(in)    :: dim        !space dimension
integer, intent(inout) :: lnn        !local number of nodes
integer, intent(in)    :: lnv        !local number of vertices
integer, intent(in)    :: lne        !local number of edges
integer, intent(in)    :: lnf        !local number of faces
integer, allocatable   :: nn(:,:)    !nodes index array
integer, allocatable   :: mm(:,:)    !vertices index array
integer :: k, j, pos
integer, allocatable :: tmp(:)

select case(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))
case('triangle')
  !insert a DOF per edge
  nparts = 0
  call alloc(tmp, size(edge_tria,1))
  do k = 1, nel
    do j = 1, lne
      tmp = sort(mm(edge_tria(:,j),k))
      call insert_sorted(1, part, tmp, used=nparts, fit=[.false.,.true.])
    end do
  end do
  call reduce(part, nparts, size(edge_tria,1))
  !add edge numeration, after vertex numeration
  lnn = lnv + lne !new nodes are vertices + edges
  call alloc(nn, lnn, nel)
  do k = 1, nel
    do j = 1, lnv !vertices
      nn(j,k) = mm(j,k)
    end do
    do j = 1, lnn-lnv !edges
      tmp = sort(mm(edge_tria(:,j),k))
      pos = find_sorted(1, part, tmp, nparts)
      if (pos > 0) nn(j+lnv,k) = nver + pos
    end do
  end do
  nnod = nver + nparts
  if (minval(nn) == 0) call error('(module_transform/lagr2p2) Some global numeration is null.')
  !re-write FE constants
  print'(a,i9)','New local number of nodes: ', lnn
  print'(a,i9)','New global number of nodes:', nnod
  call dealloc(tmp)
  call dealloc(part)

case('tetra')
  !insert a DOF per edge
  nparts = 0
  call alloc(tmp, size(edge_tetra,1))
  do k = 1, nel
    do j = 1, lne
      tmp = sort(mm(edge_tetra(:,j),k))
      call insert_sorted(1, part, tmp, used=nparts, fit=[.false.,.true.])
    end do
  end do
  call reduce(part, nparts, size(edge_tetra,1))
  !add edge numeration, after vertex numeration
  lnn = lnv + lne !new nodes are vertices + edges
  call alloc(nn, lnn, nel)
  do k = 1, nel
    do j = 1, lnv !vertices
      nn(j,k) = mm(j,k)
    end do
    do j = 1, lnn-lnv !edges
      tmp = sort(mm(edge_tetra(:,j),k))
      pos = find_sorted(1, part, tmp, nparts)
      if (pos > 0) nn(j+lnv,k) = nver + pos
    end do
  end do
  nnod = nver + nparts
  if (minval(nn) == 0) call error('(module_transform/lagr2p2) Some global numeration is null.')
  !re-write FE constants
  print'(a,i9)','New local number of nodes: ', lnn
  print'(a,i9)','New global number of nodes:', nnod
  call dealloc(tmp)
  call dealloc(part)

case ('line2', 'triangle2', 'tetra2')
  call info('(module_transform/lagr2p2) Mesh is already Lagrange P2; conversion not done.')
case default
  call error('(module_transform/lagr2p2) FE type not implemented: '//trim(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf)))
end select
call bandwidth(nel, lnn, nn, 'New Maximum bandwidth:     ')

end subroutine

!-----------------------------------------------------------------------
! lagr2rt: transforms a Lagrange FE mesh into a Lagrange Raviart-Thomas (face) FE mesh
!-----------------------------------------------------------------------
subroutine lagr2rt(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
integer, intent(in)    :: nel        !global number of elements
integer, intent(inout) :: nnod       !global number of nodes
integer, intent(inout) :: nver       !global number of vertices
integer, intent(in)    :: dim        !space dimension
integer, intent(inout) :: lnn        !local number of nodes
integer, intent(in)    :: lnv        !local number of vertices
integer, intent(in)    :: lne        !local number of edges
integer, intent(in)    :: lnf        !local number of faces
integer, allocatable   :: nn(:,:)    !nodes index array
integer, allocatable   :: mm(:,:)    !vertices index array
integer :: k, j, pos
integer, allocatable :: tmp(:)

select case(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))
case('triangle', 'triangle2')
  !insert a DOF per edge
  nparts = 0
  call alloc(tmp, size(edge_tria,1))
  do k = 1, nel
    do j = 1, lne
      tmp = sort(mm(edge_tria(:,j),k))
      call insert_sorted(1, part, tmp, used=nparts)
    end do
  end do
  call reduce(part, nparts, size(edge_tria,1))
  !add edge numeration, remove vertex numeration
  lnn = lne !new nodes are edges
  call alloc(nn, lnn, nel)
  do k = 1, nel
    do j = 1, lnn !edges
      tmp = sort(mm(edge_tria(:,j),k))
      pos = find_sorted(1, part, tmp, nparts)
      if (pos > 0) nn(j,k) = pos
    end do
  end do
  nnod = nparts
  if (minval(nn) == 0) call error('(module_transform/lagr2rt) Some global numeration is null.')
  !re-write FE constants
  print'(a,i9)','New local number of nodes: ', lnn
  print'(a,i9)','New global number of nodes:', nnod
  call dealloc(tmp)
  call dealloc(part)

case('tetra', 'tetra2')
  !insert a DOF per face
  nparts = 0
  call alloc(tmp, size(face_tetra,1))
  do k = 1, nel
    do j = 1, lnf
      tmp = sort(mm(face_tetra(:,j),k))
      call insert_sorted(1, part, tmp, used=nparts)
    end do
  end do
  call reduce(part, nparts, size(face_tetra,1))
  !add face numeration, remove vertex numeration
  lnn = lnf !new nodes are faces
  call alloc(nn, lnn, nel)
  do k = 1, nel
    do j = 1, lnn !faces
      tmp = sort(mm(face_tetra(:,j),k))
      pos = find_sorted(1, part, tmp, nparts)
      if (pos > 0) nn(j,k) = pos
    end do
  end do
  nnod = nparts
  if (minval(nn) == 0) call error('(module_transform/lagr2rt) Some global numeration is null.')
  !re-write FE constants
  print'(a,i9)','New local number of nodes: ', lnn
  print'(a,i9)','New global number of nodes:', nnod
  call dealloc(tmp)
  call dealloc(part)
case('tria-edge', 'tetra-face')
  call info('(module_transform/lagr2rt) Mesh is already Raviart-Thomas (face); conversion not done.')
case default
  call error('(module_transform/lagr2rt) FE type not implemented: '//trim(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf)))
end select
call bandwidth(nel, lnn, nn, 'New Maximum bandwidth:     ')
end subroutine

!-----------------------------------------------------------------------
! lagr2nd: transforms a Lagrange FE mesh into a Whitney (edge) FE mesh
!-----------------------------------------------------------------------
subroutine lagr2nd(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)
integer, intent(in)    :: nel        !global number of elements
integer, intent(inout) :: nnod       !global number of nodes
integer, intent(inout) :: nver       !global number of vertices
integer, intent(in)    :: dim        !space dimension
integer, intent(inout) :: lnn        !local number of nodes
integer, intent(in)    :: lnv        !local number of vertices
integer, intent(in)    :: lne        !local number of edges
integer, intent(in)    :: lnf        !local number of faces
integer, allocatable   :: nn(:,:)    !nodes index array
integer, allocatable   :: mm(:,:)    !vertices index array
integer :: k, j, pos
integer, allocatable :: tmp(:)

select case(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))
case('triangle', 'triangle2')
  !same than Raviart-Thomas
  call lagr2rt(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm)

case('tetra', 'tetra2')
  !insert a DOF per edge
  nparts = 0
  call alloc(tmp, size(edge_tetra,1))
  do k = 1, nel
    do j = 1, lne
      tmp = sort(mm(edge_tetra(:,j),k))
      call insert_sorted(1, part, tmp, used=nparts)
    end do
  end do
  call reduce(part, nparts, size(edge_tetra,1))
  !add edge numeration, remove vertex numeration
  lnn = lne !new nodes are edges
  call alloc(nn, lnn, nel)
  do k = 1, nel
    do j = 1, lne !edges
      tmp = sort(mm(edge_tetra(:,j),k))
      pos = find_sorted(1, part, tmp, nparts)
      if (pos > 0) nn(j,k) = pos
    end do
  end do
  nnod = nparts
  if (minval(nn) == 0) call error('(module_transform/lagr2nd) Some global numeration is null.')
  !re-write FE constants
  print'(a,i9)','New local number of nodes: ', lnn
  print'(a,i9)','New global number of nodes:', nnod
  call dealloc(tmp)
  call dealloc(part)
  call bandwidth(nel, lnn, nn, 'New Maximum bandwidth:     ')
case('tria-edge', 'tetra-edge')
  call info('(module_transform/lagr2nd) Mesh is already Whitney (edge); conversion not done.')
case default
  call error('(module_transform/lagr2nd) FE type not implemented: '//trim(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf)))
end select
end subroutine

!-----------------------------------------------------------------------
! to_l1: transforms a mesh into a Lagrange P1 FE mesh
!-----------------------------------------------------------------------
subroutine to_l1(pmh)
  type(pmh_mesh), intent(inout) :: pmh
  integer :: i, j, tp

  do i=1, size(pmh%pc)
    do j=1, size(pmh%pc(i)%el)
      tp = pmh%pc(i)%el(j)%type
      if(.not. FEDB(tp)%nver_eq_nnod) then
        pmh%pc(i)%el(j)%type = check_fe(.true.,FEDB(tp)%lnv, FEDB(tp)%lnv,FEDB(tp)%lne,FEDB(tp)%lnf)
        if(allocated(pmh%pc(i)%el(j)%nn)) deallocate(pmh%pc(i)%el(j)%nn)
      endif
    enddo
    pmh%pc(i)%nnod = pmh%pc(i)%nver
  enddo

end subroutine

end module
