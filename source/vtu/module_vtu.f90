module module_vtu_fcnv
!-----------------------------------------------------------------------
! Module to manage VTU files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 15/05/2013
!
! PUBLIC PROCEDURES:
!   load_vtu: loads a VTU format file
!   read_vtu: reads a mesh with fields in a VTU format file
!   save_vtu_mfm: saves a mesh in a VTU format file from mfm (refs= nsd_*, nrc_*, nra_*,nrv_*)
!   save_vtu_pmh: saves a mesh and fields in a VTU format file from pmh
!   type_cell: give the associated name of FE
!-----------------------------------------------------------------------
use basicmod
use module_pmh_fcnv
use module_fe_database_pmh_fcnv
implicit none

!Constants

integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation
!edge_tria(i,j), vertice i de la arista j de un triangulo
integer, parameter, dimension(2,3) :: edge_tria = reshape([1,2, 2,3, 3,1], [2,3])
!edge_quad(i,j), vertice i de la arista j de un cuadrangulo
integer, parameter, dimension(2,4) :: edge_quad = reshape([1,2, 2,3, 3,4, 4,1], [2,4])
!edge_tetra(i,j), vertice i de la arista j de un tetraedro
integer, parameter, dimension(2,6) :: edge_tetra = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4], [2,6])
!face_tetra(i,j), vertice i de la cara j de un tetraedro
integer, parameter, dimension(3,4) :: face_tetra = reshape([1,3,2, 1,4,3, 1,2,4, 2,3,4], [3,4])
!EDGE_HEXA(i,j), vertex #i of edge #j of a hexahedron
integer, parameter :: edge_hexa(2,12) = reshape([1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5], [2,12])
!FACE_HEXA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: face_hexa(4,6) = reshape([1,4,3,2, 1,5,8,4, 1,2,6,5, 5,6,7,8, 2,3,7,6, 3,4,8,7], [4,6])


!Private procedures
private :: true64, save_w_field, save_vtu_mfm, save_vtu_pmh !save_w_field must change

interface  save_vtu; module procedure    save_vtu_mfm; end interface
interface  save_vtu; module procedure    save_vtu_pmh; end interface

contains

!-----------------------------------------------------------------------
! save_vtu_mfm: save mesh with references
!-----------------------------------------------------------------------
subroutine save_vtu_mfm(filename, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
character(*), intent(in) :: filename
integer,      intent(in) :: nel, nnod, nver, dim, lnv, lne, lnf, lnn
integer,      intent(in) :: nn(:,:), mm(:,:), nrv(:,:), nra(:,:), nrc(:,:), nsd(:)
real(real64), intent(in) :: z(:,:)
real(real64), allocatable :: znod(:,:)
integer, dimension(:), allocatable :: ref
real(real64), allocatable, dimension(:) :: elf, elv, eln !element-wise, vertex-wise and node-wise fields
integer :: i, j, k, l, m
real(real64) :: a2(3), a3(3), a4(3)
allocate(elf(nel))
allocate(elv(nver))
allocate(eln(nnod))

!mesh
call VTU_open(trim(filename))
select case(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf)) !use P1 to store mesh
case('triangle2')
  !create znod
  allocate(znod(size(z,1), nnod))
  do k = 1, nel
    do i = 1, lnv
      znod(:,nn(i,k)) = z(:,mm(i,k))
      znod(:,nn(i+lnv,k)) = (z(:,mm(edge_tria(1,i),k))+z(:,mm(edge_tria(2,i),k)))/2
    end do
  end do
  call VTU_write_mesh(nel, nnod, nn, znod, 'triangle2')
case('tetra2')
  !create znod
  allocate(znod(size(z,1), nnod))
  do k = 1, nel
    do i = 1, lnv
      znod(:,nn(i,k)) = z(:,mm(i,k))
    end do
    do i = 1, lnn-lnv
      znod(:,nn(i+lnv, k)) = (z(:,mm(edge_tetra(1,i),k))+z(:,mm(edge_tetra(2,i),k)))/2
    end do
  end do
  call VTU_write_mesh(nel, nnod, nn, znod, 'tetra2')
case default
  call VTU_write_mesh(nel, nver, mm, z, type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))
end select

!pointdata
call VTU_begin_pointdata()
select case(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))
case('triangle') !triangles P1
  !nrv
  call sunique(pack(nrv, nrv/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1,lnv
          elv(mm(l,k)) = max(elv(mm(l,k)), true64(nrv(l,k)==ref(j)))
      end do; end do
      call VTU_write_pointdata(elv, 'nrv_'//trim(string(ref(j))), 'scalar')
    end do
  endif
  !nra
  call sunique(pack(nra, nra/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lnv; do m = 1, 2
        elv(mm(edge_tria(m,l),k)) = max(elv(mm(edge_tria(m,l),k)), true64(nra(l,k)==ref(j)))
      end do; end do; end do
      call VTU_write_pointdata(elv, 'nra_'//trim(string(ref(j))), 'scalar')
    end do
  end if
case('triangle2') !triangles P2
  !nrv
  call sunique(pack(nrv, nrv/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      eln = 0._real64
      do k = 1,nel
        do l = 1,lnv
          eln(nn(l,k)) = max(eln(nn(l,k)), true64(nrv(l,k)==ref(j)))
        end do
        do l = 1,lnn-lnv
          eln(nn(l+lnv,k)) = max(eln(nn(l+lnv,k)), (eln(nn(edge_tria(1,l),k))+eln(nn(edge_tria(2,l),k)))/2)
        end do
      end do
      call VTU_write_pointdata(eln, 'nrv_'//trim(string(ref(j))), 'scalar')
    end do
  end if
  !nra
  call sunique(pack(nra, nra/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      eln = 0._real64
      do k = 1, nel
        do l = 1, lnv; do m = 1, 2
          eln(nn(edge_tria(m,l),k)) = max(eln(nn(edge_tria(m,l),k)), true64(nra(l,k)==ref(j)))
        end do; end do
        do l = 1, lnn-lnv
          eln(nn(l+lnv,k)) = max(eln(nn(l+lnv,k)), min(eln(nn(edge_tria(1,l),k)),eln(nn(edge_tria(2,l),k))))
        end do
      end do
      call VTU_write_pointdata(eln, 'nra_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
case('quad') !Quad P1
  !nrv
  call sunique(pack(nrv, nrv/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1,lnv
          elv(mm(l,k)) = max(elv(mm(l,k)), true64(nrv(l,k)==ref(j)))
      end do; end do
      call VTU_write_pointdata(elv, 'nrv_'//trim(string(ref(j))), 'scalar')
    end do
  endif
  !nra
  call sunique(pack(nra, nra/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lnv; do m = 1, 2
        elv(mm(edge_quad(m,l),k)) = max(elv(mm(edge_quad(m,l),k)), true64(nra(l,k)==ref(j)))
      end do; end do; end do
      call VTU_write_pointdata(elv, 'nra_'//trim(string(ref(j))), 'scalar')
    end do
  end if
case('tetra') !tetrahedra P1
  !nrv
  call sunique(pack(nrv, nrv/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lnv
        elv(mm(l,k)) = max(elv(mm(l,k)), true64(nrv(l,k)==ref(j)))
      end do; end do
      call VTU_write_pointdata(elv, 'nrv_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
! nra
  call sunique(pack(nra, nra/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lne; do m = 1, 2
        elv(mm(edge_tetra(m,l),k)) = max(elv(mm(edge_tetra(m,l),k)), true64(nra(l,k)==ref(j)))
      end do; end do; end do
      call VTU_write_pointdata(elv, 'nra_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
  !nrc
  call sunique(pack(nrc, nrc/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lnf; do m = 1, 3
        elv(mm(face_tetra(m,l),k)) = max(elv(mm(face_tetra(m,l),k)), true64(nrc(l,k)==ref(j)))
      enddo; enddo; enddo
      call VTU_write_pointdata(elv, 'nrc_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
case('tetra2') !tetrahedra P2
  !nrv
  call sunique(pack(nrv, nrv/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      eln = 0._real64
      do k = 1,nel
        do l = 1,lnv
          eln(nn(l,k)) = max(eln(nn(l,k)), true64(nrv(l,k)==ref(j)))
        end do
        do l = 1, lnn-lnv
          eln(nn(l+lnv,k)) = max(eln(nn(l+lnv,k)), (eln(nn(edge_tetra(1,l),k))+eln(nn(edge_tetra(2,l),k)))/2)
        end do
      end do
      call VTU_write_pointdata(eln, 'nrv_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
! nra
  call sunique(pack(nra, nra/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      eln = 0._real64
      do k = 1, nel
        do l = 1, lne
          do m = 1, 2
            eln(nn(edge_tetra(m,l),k)) = max(eln(nn(edge_tetra(m,l),k)), true64(nra(l,k)==ref(j)))
          end do
          eln(nn(l+lnv,k))           = max(eln(nn(l+lnv,k)),           true64(nra(l,k)==ref(j)))
        end do
      end do
      call VTU_write_pointdata(eln, 'nra_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
  !nrc
  call sunique(pack(nrc, nrc/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      eln = 0._real64
      do k = 1,nel
        do l = 1, lnf; do m = 1,3
          eln(nn(face_tetra(m,l),k)) = max(eln(nn(face_tetra(m,l),k)), true64(nrc(l,k)==ref(j)))
        enddo; enddo
        do l = 1, lnn-lnv
          eln(nn(l+lnv,k)) = max(eln(nn(l+lnv,k)), min(eln(nn(edge_tetra(1,l),k)),eln(nn(edge_tetra(2,l),k))))
        end do
      enddo
      call VTU_write_pointdata(eln, 'nrc_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
case('hexahedron') !hexahedron P1
  !nrv
  call sunique(pack(nrv, nrv/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lnv
        elv(mm(l,k)) = max(elv(mm(l,k)), true64(nrv(l,k)==ref(j)))
      end do; end do
      call VTU_write_pointdata(elv, 'nrv_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
! nra
  call sunique(pack(nra, nra/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lne; do m = 1, 2
        elv(mm(edge_hexa(m,l),k)) = max(elv(mm(edge_hexa(m,l),k)), true64(nra(l,k)==ref(j)))
      end do; end do; end do
      call VTU_write_pointdata(elv, 'nra_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
  !nrc
  call sunique(pack(nrc, nrc/=0), ref)
  if (size(ref, 1) > 0) then
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1, nel; do l = 1, lnf; do m = 1, 4
        elv(mm(face_hexa(m,l),k)) = max(elv(mm(face_hexa(m,l),k)), true64(nrc(l,k)==ref(j)))
      enddo; enddo; enddo
      call VTU_write_pointdata(elv, 'nrc_'//trim(string(ref(j))), 'scalar')
    enddo
  endif

end select
call VTU_end_pointdata()

!celldata
call VTU_begin_celldata()
!nsd
call sunique(pack(nsd, nsd/=0), ref)
if (size(ref, 1) > 0) then
  do j = 1, size(ref, 1)
    where(nsd==ref(j)); elf = 1
    else where;         elf = 0
    end where
    call VTU_write_celldata(elf, 'nsd_'//trim(string(ref(j))), 'scalar')
  enddo
endif
!det
select case(type_cell(nnod, nver, dim, lnn, lnv, lne, lnf))
case('triangle', 'triangle2') !triangles
  do k = 1, nel
    a2(1:dim) = z(:,mm(2,k))-z(:,mm(1,k))
    a3(1:dim) = z(:,mm(3,k))-z(:,mm(1,k))
    elf(k) = a2(1)*a3(2)-a2(2)*a3(1)
  enddo
  call VTU_write_celldata(elf, 'det', 'scalar')
case('tetra', 'tetra2') !tetrahedrons
  do k = 1, nel
    a2 = z(:,mm(2,k))-z(:,mm(1,k))
    a3 = z(:,mm(3,k))-z(:,mm(1,k))
    a4 = z(:,mm(4,k))-z(:,mm(1,k))
    elf(k) = a2(1)*a3(2)*a4(3) + a2(3)*a3(1)*a4(2) + a2(2)*a3(3)*a4(1) &
            -a2(3)*a3(2)*a4(1) - a2(2)*a3(1)*a4(3) - a2(1)*a3(3)*a4(2)
  enddo
  call VTU_write_celldata(elf, 'det', 'scalar')
end select
call VTU_end_celldata()
call VTU_close()

end subroutine


!-----------------------------------------------------------------------
! save_vtu_pmh(filename, pmh): write VTU file
!-----------------------------------------------------------------------
! filename:   name of a VTU file
! pmh:    PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------
subroutine save_vtu_pmh(filename, pmh, infield, outfield, padval, nparam, param, cell2point)
character(len=*),              intent(in)  :: filename
type(pmh_mesh),                intent(inout)  :: pmh ! pmh_mesh
character(len=*), allocatable, intent(in)  :: infield(:)
character(len=*), allocatable, intent(in)  :: outfield(:)
real(real64),                  intent(in)  :: padval
integer,          optional,    intent(in)  :: nparam
real(real64),     optional,    intent(out) :: param
logical,          optional,    intent(in)  :: cell2point
character(len=maxpath) :: fieldname, str
integer :: i, j, k, l, m, nnod, nel, tp, lnn, lnv, tnvpc, maxtopdim, np
logical :: all_P1
real(real64), allocatable :: znod(:,:), pdfvalmat2(:,:)
integer, allocatable :: connect(:), offset(:), celltypes(:),piece2save(:), v_ref(:), e_ref(:), f_ref(:), el_ref(:), uref(:), &
                        aux_ref(:), temp(:,:)
type cdfield
  character(len=maxpath)    :: name
  integer                   :: ncomp
  real(real64), allocatable :: val(:,:)
end type
type(cdfield), allocatable    :: cdfval(:)

if(allocated(outfield)) then
  if(allocated(infield)) then
    if(size(outfield,1) /= size(infield,1)) call error('Number of input/output field names must agree.')
  elseif(size(outfield,1) /= get_piece_num_fields(pmh%pc(i))) then
    call error('Number of field names must agree with the number of fields.')
  endif
endif

!check piece(s) to be saved
if (is_arg('-p')) then !save a single piece, indicated after -p
  str = get_post_arg('-p')
  call alloc(piece2save, 1)
  piece2save(1) = int(str)
else !save all pieces
  call alloc(piece2save, size(pmh%pc,1))
  piece2save = [(i, i=1, size(pmh%pc,1))]
end if
np = 1
if(present(nparam)) np = nparam
call VTU_open(trim(filename))
! Loop in pieces
call info('Number of pieces: '//trim(string(size(pmh%pc,1))))
do i=1,size(piece2save,1)
  call info('  Piece: '//trim(string(piece2save(i))))
  nel = 0   ! Number of VTK elements
  tnvpc = 0 ! Total number of vertex per cell
  maxtopdim = 0
  associate(pc => pmh%pc(piece2save(i)))
    call build_node_coordinates(pc, i, all_P1, znod)
    if(allocated(v_ref)) deallocate(v_ref)
    if(.not. all_P1)then
      allocate(v_ref(size(znod,2)))
    else
      allocate(v_ref(pc%nver))
    endif
    v_ref = 0
    ! Calc max topological dimension
    do j=1,size(pc%el,1)
      maxtopdim = max(maxtopdim,FEDB(pc%el(j)%type)%tdim)
    enddo
    ! Loop in element groups
    do j = 1, size(pc%el,1)
      tp = pc%el(j)%type
      call info('    Element type: '//trim(FEDB(tp)%desc))
      if(.not. FEDB(tp)%nver_eq_nnod) then
        lnn = FEDB(tp)%lnn
        ! Build vtk conectivity array
        if (tp == check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron P2 needs node reordering
          ! VTU HexaP2 nodes: [1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16]
          call set(connect, &
          & pack(reshape(pc%el(j)%nn([1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16],:)-1,&
          & [1,pc%el(j)%nel*lnn]),.true.), &
          & [(k,k=tnvpc+1,tnvpc+pc%el(j)%nel*lnn)], fit=.false.)
        else
          ! Build vtk conectivity array
          call set(connect, pack(reshape(pc%el(j)%nn-1, [1,pc%el(j)%nel*lnn]),.true.), &
          &[(k,k=tnvpc+1,tnvpc+pc%el(j)%nel*lnn)], fit=.false.)
        endif
        ! Build vtk offset array
        call set(offset, [(tnvpc+k*lnn, k=1,pc%el(j)%nel)], [(k, k = nel+1,nel+pc%el(j)%nel)], fit=.false.)
        tnvpc = tnvpc + pc%el(j)%nel*lnn
      else
        lnv = FEDB(tp)%lnv
        ! Build vtk conectivity array
        if(allocated(temp)) deallocate(temp)
        !allocate(temp, source=reshape(pc%el(j)%mm-1, [1,pc%el(j)%nel*lnv]))
        allocate(temp(1,pc%el(j)%nel*lnv))
        temp(:,:) = reshape(pc%el(j)%mm-1, [1, pc%el(j)%nel*lnv])
        call set(connect, temp(1,:), [(k, k = tnvpc+1,tnvpc+pc%el(j)%nel*lnv)], fit=.false.)
        deallocate(temp)
        ! Build vtk offset array
        call set(offset, [(tnvpc+k*lnv, k = 1,pc%el(j)%nel)], [(k,k=nel+1,nel+pc%el(j)%nel)], fit=.false.)
        tnvpc = tnvpc + pc%el(j)%nel*lnv
      endif
      ! Build vtk cell types array
      call set(celltypes, [(pmh2vtkcelltype(tp), k=nel,nel+pc%el(j)%nel)] , [(k, k = nel+1,nel+1+pc%el(j)%nel)], fit=.false.)
      if(FEDB(tp)%tdim == maxtopdim) then !subdomain references
        call set(el_ref, [(pc%el(j)%ref(k), k=1,pc%el(j)%nel)] , [(k, k = nel+1,nel+1+pc%el(j)%nel)], fit=.false.)
      elseif(FEDB(tp)%tdim == 2) then     !face references
        call set(f_ref, [(pc%el(j)%ref(k), k=1,pc%el(j)%nel)] , [(k, k = nel+1,nel+1+pc%el(j)%nel)], fit=.false.)
      elseif(FEDB(tp)%tdim == 1) then     !edge references
        call set(e_ref, [(pc%el(j)%ref(k), k=1,pc%el(j)%nel)] , [(k, k = nel+1,nel+1+pc%el(j)%nel)], fit=.false.)
      elseif(FEDB(tp)%tdim == 0) then     !vertex references
        do k=1,pc%el(j)%nel
          v_ref(pc%el(j)%mm(:,k)) = pc%el(j)%ref(k)
        enddo
        !call set(v_ref, [(pc%el(j)%ref(k), k=1,pc%el(j)%nel)], [(k,k=nel+1,nel+1+pc%el(j)%nel)], fit=.false.) !vertex as celldata
      endif
      if(.not. allocated(cdfval) .and. allocated(pc%el(j)%fi)) allocate(cdfval(size(pc%el(j)%fi,1)))
      if(allocated(pc%el(j)%fi)) then
        do l=1, size(pc%el(j)%fi,1)
          cdfval(l)%name = trim(pc%el(j)%fi(l)%name)
          cdfval(l)%ncomp = size(pc%el(j)%fi(l)%val,1)
          if(.not. allocated(cdfval(l)%val)) then
            allocate(cdfval(l)%val(cdfval(l)%ncomp, nel+pc%el(j)%nel))
            cdfval(l)%val =  padval
          endif
          do m=1, size(pc%el(j)%fi(l)%val,1)
            call alloc(pdfvalmat2, 1, pc%el(j)%nel)
            pdfvalmat2 = reshape([(pc%el(j)%fi(l)%val(m,k,np), k=1,pc%el(j)%nel)], [1, pc%el(j)%nel])
            call set(cdfval(l)%val, pdfvalmat2, [m], [(k,k=nel+1,nel+pc%el(j)%nel)])
          enddo
          if(present(param)) param = pc%el(j)%fi(l)%param(np)
        enddo
      endif
      nel = nel + pc%el(j)%nel
    enddo
  
    if(allocated(connect))   call reduce(connect,   tnvpc)
    if(allocated(offset))    call reduce(offset,    nel)
    if(allocated(celltypes)) call reduce(celltypes, nel)
    if(pc%dim < 3) then
      if(.not. all_P1) then
        nnod = size(znod,2)
        if (vtk_geo_xml(nnod, nel, pack(znod(1,:),.true.), pack(znod(2,:),.true.), [(0._real64,l=1,nnod)]) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      else
        nnod = pc%nver
        if (vtk_geo_xml(pc%nver, nel, &
        &  pack(pc%z(1,:),.true.), pack(pc%z(2,:),.true.), [(0._real64,l=1,pc%nver)]) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      endif
    else
      if(.not. all_P1) then
        nnod = size(znod,2)
        if (vtk_geo_xml(nnod, nel, pack(znod(1,:),.true.), pack(znod(2,:),.true.), pack(znod(3,:),.true.)) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      else
        nnod = pc%nver
        if (vtk_geo_xml(pc%nver, nel, &
        &  pack(pc%z(1,:),.true.), pack(pc%z(2,:),.true.), pack(pc%z(3,:),.true.)) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      endif
    endif
  
    if(vtk_con_xml(nel,connect,offset,int(celltypes,I1P)) /= 0) stop 'writeVTU: vtk_con_xml error 1'
    call VTU_begin_pointdata()
    ! Array used to save references in this way: nsd_*
    if(allocated(aux_ref)) deallocate(aux_ref)
    allocate(aux_ref(nnod))
    if(allocated(v_ref) .and. size(pack(v_ref,v_ref/=0))>0) then
      call info('    Writing vertex references')
      if(size(v_ref,1)<nnod) call add(v_ref,[(0,k=size(v_ref,1),nnod)],[(k,k=size(v_ref,1),nnod)],fit=.true.)
      call reduce(v_ref, nnod)
      if(vtk_var_xml(nnod, 'vertex_ref', v_ref) /= 0) stop
      if(allocated(uref)) deallocate(uref)
      call sunique(pack(v_ref,v_ref /= 0), uref)
       uref = unique(pack(v_ref,v_ref /= 0))
      do k=1, size(uref,1)
        aux_ref = 0
        where(v_ref==uref(k))
          aux_ref = v_ref
        end where
        if(vtk_var_xml(nnod, 'nrv_'//trim(string(uref(k))), aux_ref) /= 0) stop
      enddo
    endif
  
    if(allocated(pc%fi)) then
      do j=1, size(pc%fi,1)
        if(allocated(outfield)) then
          if(allocated(infield)) then
            do k=1,size(infield,1)
              if(trim(pc%fi(j)%name) == infield(k)) fieldname = outfield(k)
            enddo
          elseif(size(outfield,1) == get_piece_num_fields(pc)) then
            fieldname = outfield(k)
          else
            call error('Number of field names must agree with the number of fields.')
          endif
        else
         fieldname = trim(pc%fi(j)%name)
        endif
        call info('    Writing node field: '//trim(fieldname))
        if(size(pc%fi(j)%val,1) == 1) then
          if (vtk_var_xml(nnod, trim(fieldname), pc%fi(j)%val(1,1:nnod,np)) /= 0) call error('Writing '//trim(fieldname))
        elseif(size(pc%fi(j)%val,1) == 2) then
          if (vtk_var_xml(nnod, trim(fieldname), &
          pc%fi(j)%val(1,1:nnod,np), pc%fi(j)%val(2,1:nnod,np), [(real(0,R8P),i=1,nnod)]) /= 0) &
          call error('Writing '//trim(fieldname))
        elseif(size(pc%fi(j)%val,1) == 3) then
          if (vtk_var_xml(nnod, trim(fieldname), &
          real(pc%fi(j)%val(1,1:nnod,np),R8P), real(pc%fi(j)%val(2,1:nnod,np),R8P), &
          real(pc%fi(j)%val(3,1:nnod,np),R8P)) /= 0) call error('Writing '//fieldname)
        else
          if( vtk_var_xml(nnod,size(pc%fi(j)%val,1),trim(fieldname), &
          reshape(transpose(pc%fi(j)%val(:,1:nnod,np)),[nnod, size(pc%fi(j)%val,1)]) ) /= 0) &
          call error('Writing '//fieldname)
        endif
        if(present(param)) param = pc%fi(j)%param(np)
      enddo
    endif
    call VTU_end_pointdata()
    call VTU_begin_celldata()
    ! Array used to save references in this way: nsd_*
    if(allocated(aux_ref)) deallocate(aux_ref)
    allocate(aux_ref(nel))
  
    if(allocated(el_ref) .and. size(pack(el_ref,el_ref/=0))>0) then
      call info('    Writing subdomain references')
      if(size(el_ref,1)<nel) call add(el_ref,[(0,k=size(el_ref,1),nel)],[(k,k=size(el_ref,1),nel)],fit=.true.)
      call reduce(el_ref, nel)
      if(vtk_var_xml(nel, 'element_ref', el_ref) /= 0) stop
      ! Build nsd references
      if(allocated(uref)) deallocate(uref)
      call sunique(pack(el_ref,el_ref /= 0), uref)
      do k=1, size(uref,1)
        aux_ref = 0
        where(el_ref==uref(k))
          aux_ref = el_ref
        end where
        if(vtk_var_xml(nel, 'nsd_'//trim(string(uref(k))), aux_ref) /= 0) stop
      enddo
    endif
    if (allocated(f_ref)) then
      call info('    Writing face references')
      if (size(f_ref,1) < nel) call add(f_ref, [(0, k = size(f_ref,1), nel)], [(k, k = size(f_ref,1),nel)], fit=.true.)
      call reduce(f_ref, nel)
      if(vtk_var_xml(nel, 'face_ref', f_ref) /= 0) stop
      ! Build nrc references
      if (allocated(uref)) deallocate(uref)
      call sunique(pack(f_ref,f_ref /= 0), uref)
      do k = 1, size(uref, 1)
        aux_ref = 0
        where(f_ref == uref(k))
          aux_ref = f_ref
        end where
        if (vtk_var_xml(nel, 'nrc_'//trim(string(uref(k))), aux_ref) /= 0) stop
      enddo
    endif
    if (allocated(e_ref)) then
      call info('    Writing edge references')
      if(size(e_ref,1)<nel) call add(e_ref,[(0,k=size(e_ref,1),nel)],[(k,k=size(e_ref,1),nel)],fit=.true.)
      call reduce(e_ref, nel)
      if(vtk_var_xml(nel, 'edge_ref', e_ref) /= 0) stop
      ! Build nra references
      if(allocated(uref)) deallocate(uref)
      call sunique(pack(e_ref,e_ref /= 0),uref)
       uref = unique(pack(e_ref,e_ref /= 0))
      do k=1, size(uref,1)
        aux_ref = 0
        where(e_ref==uref(k))
          aux_ref = e_ref
        end where
        if(vtk_var_xml(nel, 'nra_'//trim(string(uref(k))), aux_ref) /= 0) stop
      enddo
    endif
  
    if(allocated(cdfval)) then
      do j=1, size(cdfval,1)
        if(allocated(outfield)) then
          if(allocated(infield)) then
            do k=1,size(infield,1)
              if(trim(cdfval(j)%name) == infield(k)) fieldname = outfield(k)
            enddo
          elseif(size(outfield,1) == get_piece_num_fields(pc)) then
            fieldname = outfield(k)
          endif
        else
          fieldname = trim(cdfval(j)%name)
        endif
        call info('    Writing cell field: '//trim(fieldname))
        if(size(cdfval(j)%val,2) < nel) then
          call alloc(pdfvalmat2, size(cdfval(j)%val,1), nel-size(cdfval(j)%val,2)+1)
          pdfvalmat2 = reshape([(padval, k=size(cdfval(j)%val,2),nel*size(cdfval(j)%val,1))], &
                               [size(cdfval(j)%val,1), nel-size(cdfval(j)%val,2)+1])
          call add(cdfval(j)%val, pdfvalmat2, [(k,k=1,size(cdfval(j)%val,1))], [(k,k=size(cdfval(j)%val,2),nel)])
        end if
        call reduce(cdfval(j)%val, cdfval(j)%ncomp, nel)
       if(size(cdfval(j)%val,1) == 1) then
          if (vtk_var_xml(nel, trim(fieldname), reshape(cdfval(j)%val(1,1:nel), [nel])) /= 0) &
          call error('Writing '//trim(fieldname))
        elseif(size(cdfval(j)%val,1) == 2) then
          if (vtk_var_xml(nel, trim(fieldname), &
          reshape(cdfval(j)%val(1,1:nel),[nel]), reshape(cdfval(j)%val(2,1:nel),[nel]), [(real(0,R8P),i=1,nel)]) /= 0) &
          call error('Writing '//trim(fieldname))
        elseif(size(cdfval(j)%val,1) == 3) then
          if (vtk_var_xml(nel, trim(fieldname), &
          reshape(cdfval(j)%val(1,1:nel),[nel]), reshape(cdfval(j)%val(2,1:nel),[nel]), &
          reshape(cdfval(j)%val(3,1:nel),[nel])) /= 0) call error('Writing '//trim(fieldname))
        else
          if( vtk_var_xml(nel, size(cdfval(j)%val,1), trim(fieldname),&
          reshape(transpose(cdfval(j)%val(:,:)),[nel, size(cdfval(j)%val,1)]) ) /= 0) &
          call error('Writing '//fieldname)
          !call error('Vector field with '//trim(string(size(cdfval(j)%val,1)))//' components not supported!')
        endif
      enddo
    endif
    call VTU_end_celldata()
    if(vtk_geo_xml() /= 0) stop 'writeVTU: vtk_geo_xml_closep error 1'
  end associate
enddo
call VTU_close(geo=.false.)
end subroutine

!-----------------------------------------------------------------------
! save_w_field: save mesh with field
!-----------------------------------------------------------------------
subroutine save_w_field(filename, nel, nnod, nver, lnn, dim, lnv, lne, lnf, &
nn, mm, z, fichfield, v)
character(len=*), intent(in)              :: filename, fichfield
integer,          intent(in)              :: nel, nnod, nver, dim, lnv, lne, lnf, lnn
integer,       dimension(:,:) :: nn, mm
real(real64),  dimension(:,:) :: z(:,:), v(:)
real(real64), allocatable :: znod(:,:)
integer :: i, k
character(80) :: tc

tc = type_cell(nnod, nver, dim, lnn, lnv, lne, lnf)
call VTU_open(trim(filename))
!triangles P2 with field P2
!----------------------------------------------------------
if (tc == 'triangle2' .and. mod(size(v,1), nnod) == 0) then
  !create znod
  allocate(znod(size(z,1), nnod))
  do k = 1, nel
    do i = 1, 3
      znod(:,nn(i,k)) = z(:,mm(i,k))
    enddo
    znod(:,nn(4,k)) = (z(:,mm(1,k))+z(:,mm(2,k)))/2
    znod(:,nn(5,k)) = (z(:,mm(2,k))+z(:,mm(3,k)))/2
    znod(:,nn(6,k)) = (z(:,mm(3,k))+z(:,mm(1,k)))/2
  enddo
  call VTU_write_mesh(nel, nnod, nn, znod, tc)
  call VTU_begin_pointdata()
  select case(size(v,1)/nnod) !v components
  case(1) !scalar field
    call VTU_write_pointdata(v, trim(fichfield), 'scalar')
  case(2) !2D field
    if (size(z,1) /= 2) call error('Incorrect field dimension.')
    call VTU_write_pointdata(v, trim(fichfield), 'vector')
  case(3) !3D field
    if (size(z,1) /= 3) call error('Incorrect field dimension.')
    call VTU_write_pointdata(v, trim(fichfield), 'vector')
  end select
  call VTU_end_pointdata()
!triangles with field P1
!----------------------------------------------------------
elseif ((tc == 'triangle2' .or. tc == 'triangle') .and. mod(size(v,1), nver) == 0) then
  tc = 'triangle'
  call VTU_write_mesh(nel, nver, mm, z, tc)
  call VTU_begin_pointdata()
  select case(size(v,1)/nver) !v components
  case(1) !scalar field
    call VTU_write_pointdata(v, trim(fichfield), 'scalar')
  case(2) !2D field
    if (size(z,1) /= 2) call error('Incorrect field dimension.')
    call VTU_write_pointdata(v, trim(fichfield), 'vector')
  case(3) !3D field
    if (size(z,1) /= 3) call error('Incorrect field dimension.')
    call VTU_write_pointdata(v, trim(fichfield), 'vector')
  end select
  call VTU_end_pointdata()
!tetrahedra P2 with field P2
!----------------------------------------------------------
elseif (tc == 'tetra2' .and. mod(size(v,1), nnod) == 0) then
  !create znod
  allocate(znod(size(z,1), nnod))
  do k = 1, nel
    do i = 1, 4
      znod(:,nn(i,k)) = z(:,mm(i,k))
    enddo
    znod(:,nn(5, k)) = (z(:,mm(1,k))+z(:,mm(2,k)))/2
    znod(:,nn(6, k)) = (z(:,mm(2,k))+z(:,mm(3,k)))/2
    znod(:,nn(7, k)) = (z(:,mm(3,k))+z(:,mm(1,k)))/2
    znod(:,nn(8, k)) = (z(:,mm(1,k))+z(:,mm(4,k)))/2
    znod(:,nn(9, k)) = (z(:,mm(2,k))+z(:,mm(4,k)))/2
    znod(:,nn(10,k)) = (z(:,mm(3,k))+z(:,mm(4,k)))/2
  enddo
  call VTU_write_mesh(nel, nnod, nn, znod, tc)
  call VTU_begin_pointdata()
  select case(size(v,1)/nnod)
  case(1) !scalar field
    call VTU_write_pointdata(v, trim(fichfield), 'scalar')
  case(3) !3D field
    call VTU_write_pointdata(v, trim(fichfield), 'vector')
  end select
  call VTU_end_pointdata()
!tetrahedra with field P1
!----------------------------------------------------------
elseif ((tc == 'tetra2' .or. tc == 'tetra') .and. mod(size(v,1), nver) == 0) then
  tc = 'tetra'
  call VTU_write_mesh(nel, nver, mm, z, tc)
  call VTU_begin_pointdata()
  select case(size(v,1)/nver)
  case(1) !scalar field
    call VTU_write_pointdata(v, trim(fichfield), 'scalar')
  case(3) !3D field
    call VTU_write_pointdata(v, trim(fichfield), 'vector')
  end select
  call VTU_end_pointdata()
!celldata
!----------------------------------------------------------
elseif (size(v,1) == nel) then
  tc = 'tetra'
  call VTU_write_mesh(nel, nver, mm, z, tc)
  call VTU_begin_celldata()
  call VTU_write_celldata(v, trim(fichfield), 'scalar')
  call VTU_end_celldata()
else
  call error('Element not supported.')
end if
call VTU_close()

end subroutine

!-----------------------------------------------------------------------
! type_cell: give the associated name of FE
!-----------------------------------------------------------------------
function type_cell(nnod, nver, DIM, LNN, LNV, LNE, LNF) result(res)
integer, intent(in) :: nnod  !total number of nodes
integer, intent(in) :: nver  !total number of vertices
integer, intent(in) :: DIM   !space dimension
integer, intent(in) :: LNN   !local number of nodes
integer, intent(in) :: LNV   !local number of vertices
integer, intent(in) :: LNE   !local number of edges
integer, intent(in) :: LNF   !local number of faces
character(10)       :: res   !FE name

if     (             LNV==2 .and. LNN== 2 .and. LNE<= 1 .and. nver == nnod) then
  res = 'line'
elseif (             LNV==2 .and. LNN== 3 .and. LNE<= 1) then
  res = 'line2'
elseif (             LNV==3 .and. LNN== 3 .and. LNE== 3 .and. nver == nnod) then
  res = 'triangle'
elseif (             LNV==3 .and. LNN== 6 .and. LNE== 3) then
  res = 'triangle2'
elseif (             LNV==3 .and. LNN== 3 .and. LNE== 3 .and. nver /= nnod) then
  res = 'tria-edge'
elseif (             LNV==4 .and. LNN== 4 .and. LNE== 4 .and. nver == nnod) then
  res = 'quad'
elseif (DIM==3 .and. LNV==4 .and. LNN== 4 .and. LNE== 6 .and. LNF==4 .and. nver == nnod) then
  res = 'tetra'
elseif (DIM==3 .and. LNV==4 .and. LNN==10 .and. LNE== 6 .and. LNF==4) then
  res = 'tetra2'
elseif (DIM==3 .and. LNV==4 .and. LNN== 6 .and. LNE== 6 .and. LNF==4) then
  res = 'tetra-edge'
elseif (DIM==3 .and. LNV==4 .and. LNN==20 .and. LNE== 6 .and. LNF==4) then
  res = 'tetr-edge2'
elseif (DIM==3 .and. LNV==4 .and. LNN== 4 .and. LNE== 6 .and. LNF==4 .and. nver /= nnod) then
  res = 'tetra-face'
elseif (DIM==3 .and. LNV==8 .and. LNN== 8 .and. LNE==12 .and. LNF==6) then
  res = 'hexahedron'
elseif (DIM==3 .and. LNV==6 .and. LNN== 6 .and. LNE== 9 .and. LNF==5) then
  res = 'wedge'
elseif (DIM==3 .and. LNV==5 .and. LNN== 5 .and. LNE== 8 .and. LNF==5) then
  res = 'pyramid'
else
  call error('(module_vtu/type_cell), cell not recognized: DIM = '//trim(string(DIM))//', LNV = '//trim(string(LNV))//', LNN = '//&
  &trim(string(LNN))//', LNE = '//trim(string(LNE))//', LNF = '//trim(string(LNF))//', nver==nnod: '//trim(string(nver == nnod)))
end if
end function

!-----------------------------------------------------------------------
! true64: 1. when true, 0. when false
!-----------------------------------------------------------------------
pure real(8) function true64(l)
logical, intent(in) :: l
true64 = 0._8
if (l) true64 = 1._8
end function

subroutine load_vtu(filename, pmh, infieldname, nparam, param)
  character(len=*),             intent(in) :: filename
  type(pmh_mesh),            intent(inout) :: pmh
  character(len=*), allocatable,intent(in) :: infieldname(:)
  integer, optional,            intent(in) :: nparam
  real(real64), optional,       intent(in) :: param

  ! Inital settings
  !call report_option('level', 'stdout')

  ! Open VTU file and reads the mesh
  call info('Reading VTU file ...')
  call read_vtu(filename, pmh, infieldname, nparam, param)

  call build_vertices(pmh)

end subroutine


!-----------------------------------------------------------------------
! read_VTU(filename, pmh): read VTU file
!-----------------------------------------------------------------------
! filename:   VTU file name
! pmh:    PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------
subroutine read_vtu(filename, pmh, fieldnames, nparam, param)


  type cdfield
    real(R8P), allocatable    :: val(:,:)
  end type

  character(len=*),              intent(in) :: filename
  type(pmh_mesh),             intent(inout) :: pmh ! pmh_mesh
  character(len=*), allocatable, intent(in) :: fieldnames(:)
  integer, optional,             intent(in) :: nparam
  real(real64), optional,        intent(in) :: param

  character(len=MAXPATH), allocatable       :: pdfnames(:), cdfnames(:)
  integer(I4P)              :: ip, nn, nc, nnref, ncref, ncomp
  real(R8P), allocatable    :: X(:), Y(:), Z(:), pdfval(:), pdfvalmat3(:,:,:)
  real(R8P), allocatable    :: temp(:),temp2(:,:,:)
  integer, allocatable      :: v_ref(:), aux_ref(:), c_refs(:)
  integer(I4P), allocatable :: connect(:), offset(:)
  integer(I1P), allocatable :: cell_type(:),uct(:)
  integer                   :: i,j,k, l,lnn, npdf, ncdf, np, nel
  type(cdfield), allocatable:: cdfval(:)
  real(R8P)                 :: p
  logical                   :: file_exists, ffound
  logical                   :: tf(2)=[.true.,.false.],tt(2)=[.true.,.true.]
  type(pmh_mesh)            :: auxpmh
  logical                   :: ttt(3) = [.true.,.true.,.true.]


  ! Reads the VTU file header and allocates the number of pieces
  inquire(file=trim(filename), exist=file_exists)

  if(.not. file_exists) call error('Input file '//trim(filename)//' not found!')

  if (vtk_ini_xml_read('Binary',filename,'UnstructuredGrid', ip)/=0) stop 'Error'
  call info('Number of pieces: '//trim(string(ip)))

  np = 1
  if(present(nparam)) np = nparam
  p = 0._real64
  if(present(param)) p = param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if(allocated(pmh%pc)) deallocate(pmh%pc)
!  allocate(pmh%pc(ip))

  ! Memory allocation: Pieces
  if(.not. allocated(pmh%pc)) then
    allocate(pmh%pc(ip))
  elseif(size(pmh%pc,1) < ip) then
    call move_alloc(from=pmh%pc, to=auxpmh%pc)
    if(allocated(pmh%pc)) deallocate(pmh%pc)
    allocate(pmh%pc(ip))
    pmh%pc(1:size(auxpmh%pc,1)) = auxpmh%pc(1:size(auxpmh%pc,1))
    deallocate(auxpmh%pc)
  endif

  if(.not. allocated(auxpmh%pc)) allocate(auxpmh%pc(1))
  do i=1, ip
    call info('  Piece '//trim(string(i))//':')
    call info('    Reading coordinates ...')
    if (vtk_geo_xml_read(nn, nc, X, Y, Z, i)/=0) stop 'Error'
    pmh%pc(i)%nnod = nn
    pmh%pc(i)%dim = 3
    if(allocated(pmh%pc(i)%z)) deallocate(pmh%pc(i)%z)
    allocate(pmh%pc(i)%z(3,nn))

    ! Save coordinates in PMH structure
    call set(1, pmh%pc(i)%z, X, 1, fit=tt)
    if(allocated(X)) deallocate(X)
    call set(1, pmh%pc(i)%z, Y, 2, fit=tt)
    if(allocated(Y)) deallocate(Y)
    call set(1, pmh%pc(i)%z, Z, 3, fit=tt)
    if(allocated(Z)) deallocate(Z)


    call info('    Reading conectivities ... ')
    if (vtk_con_xml_read(nc,connect,offset,cell_type,i)/=0) stop 'Error'

    ! Total array of cell references (celldata)
    if(allocated(c_refs)) deallocate(c_refs)
    allocate(c_refs(nc))
    c_refs = 0


    call info('    Reading references ... ')
    ! VTK+ vertex references (pointdata)
    if(vtk_var_xml_read('Node',nnref,ncomp,'vertex_ref',v_ref,i) == 0 ) then
      if((nnref/=nn .or. ncomp /= 1) .and. allocated(v_ref)) deallocate(v_ref)
    endif
    ! VTK+ edge references (celldata)
    if(vtk_var_xml_read('Cell',ncref,ncomp,'edge_ref',aux_ref,i) == 0) then
      if((ncref/=nc .or. ncomp/= 1) .and. allocated(aux_ref)) then
        deallocate(aux_ref)
      elseif(allocated(aux_ref)) then
         call add(c_refs,aux_ref,[(j,j=1,nc)],fit=.true.)
         deallocate(aux_ref)
      endif
    endif
    ! VTK+ face references (celldata)
    if(vtk_var_xml_read('Cell',ncref,ncomp,'face_ref',aux_ref,i) == 0) then
      if((ncref/=nc .or. ncomp/=1) .and. allocated(aux_ref)) then
        deallocate(aux_ref)
      elseif(allocated(aux_ref)) then
        call add(c_refs,aux_ref,[(j,j=1,nc)],fit=.true.)
        deallocate(aux_ref)
      endif
    endif
    ! VTK+ element references (celldata)
    if(vtk_var_xml_read('Cell',ncref,ncomp,'element_ref',aux_ref,i) == 0) then
      if((ncref/=nc .or. ncomp /= 1) .and. allocated(aux_ref)) then
        deallocate(aux_ref)
      elseif(allocated(aux_ref)) then
        call add(c_refs,aux_ref,[(j,j=1,nc)],fit=.true.)
        deallocate(aux_ref)
      endif
    endif

    ! Reads and store pointdata fields
    if (vtk_var_xml_list(pdfnames,i,'Node') == 0 ) then
      if (allocated(pdfnames)) then
        npdf = size(pdfnames,1) ! Number of pointdata fields
        do j=1, size(pdfnames,1)
          if (trim(lcase(pdfnames(j))) == 'vertex_ref' .or. index(lcase(pdfnames(j)),'nrv_') /= 0) then
            npdf = npdf - 1
          elseif (allocated(fieldnames)) then !if -in option is present, discard not selected fields
            ffound = .false.
            do k=1, size(fieldnames,1)
              if (trim(adjustl(pdfnames(j))) == trim(adjustl(fieldnames(k)))) then
                ffound = .true.
                exit
              end if  
            end do
            if(.not. ffound) npdf = npdf - 1
          end if
        end do
      else
        npdf = 0
      end if

      ! Memory allocation: Point fields
      if(.not. allocated(pmh%pc(i)%fi)) then
        allocate(pmh%pc(i)%fi(npdf))
      elseif(size(pmh%pc(i)%fi,1) < npdf) then
        call move_alloc(from=pmh%pc(i)%fi, to=auxpmh%pc(1)%fi)
        if(allocated(pmh%pc(i)%fi)) deallocate(pmh%pc(i)%fi)
        allocate(pmh%pc(i)%fi(npdf))
        pmh%pc(i)%fi(1:size(auxpmh%pc(1)%fi,1)) = auxpmh%pc(1)%fi(:)
        if(allocated(auxpmh%pc(1)%fi)) deallocate(auxpmh%pc(1)%fi)
      endif

      npdf = 0
      if(allocated(fieldnames) .and. allocated(pdfnames)) then
        do k=1, size(fieldnames,1)
          do j=1, size(pdfnames,1)
            if(trim(adjustl(pdfnames(j))) == trim(adjustl(fieldnames(k)))) then
              if(vtk_var_xml_read('Node', nnref,ncomp, trim(pdfnames(j)), pdfval, i) == 0 ) then
                npdf = npdf +1
                if(nnref == nn) then
                  call info('    Reading node field: '//trim(pdfnames(j)))
                  pmh%pc(i)%fi(npdf)%name = trim(pdfnames(j))
                  call  set(pmh%pc(i)%fi(npdf)%param, p, np, fit=.true.)
                  call alloc(pdfvalmat3, ncomp, nnref, 1)
                  pdfvalmat3 = reshape(pdfval, [ncomp, nnref, 1])
                  call set(pmh%pc(i)%fi(npdf)%val, pdfvalmat3, [(k,k=1,ncomp)], [(k,k=1,nnref)], [np], fit=ttt)
                  call reduce(pmh%pc(i)%fi(npdf)%val, ncomp, nnref, np)
                else
                  call error('Wrong number of node values')
                endif
              else
                call error('Reading VTU file')
              endif
            endif
          enddo
        enddo
      elseif(allocated(pdfnames)) then
        do j=1, size(pdfnames,1)
          if(lcase(pdfnames(j)) == 'vertex_ref' .or. &
            & index(lcase(pdfnames(j)),'nrv_') /= 0) cycle
          ! if -in option is present, discard not selected fields
          if(vtk_var_xml_read('Node', nnref, ncomp, trim(pdfnames(j)), pdfval, i) == 0 ) then
            npdf = npdf +1
            if(nnref == nn) then
              call info('    Reading node field: '//trim(pdfnames(j)))
              pmh%pc(i)%fi(npdf)%name = trim(pdfnames(j))
              call set(pmh%pc(i)%fi(npdf)%param, p, np, fit=.true.)
              call alloc(pdfvalmat3, ncomp, nnref, 1)
              pdfvalmat3 = reshape(pdfval, [ncomp, nnref, 1])
              call set(pmh%pc(i)%fi(npdf)%val, pdfvalmat3, [(k,k=1,ncomp)], [(k,k=1,nnref)], [np], fit=ttt)
              call reduce(pmh%pc(i)%fi(npdf)%val, ncomp, nnref, np)
            else
              call error('Wrong number of node values')
            endif
          else
            call error('Reading VTU file')
          endif
        enddo
      endif
    endif
    if(allocated(pdfval)) deallocate(pdfval)
    if(allocated(pdfnames)) deallocate(pdfnames)

    ncdf = 0
    ! Reads celldata fields
    if(vtk_var_xml_list(cdfnames,i,'Cell') == 0 ) then
      if(allocated(cdfnames)) then
        ncdf = size(cdfnames,1) ! Number of celldata fields
        do j=1, size(cdfnames,1)
          ! Discard references
          if (trim(lcase(cdfnames(j))) == 'element_ref' .or. trim(lcase(cdfnames(j))) == 'face_ref' .or. &
              trim(lcase(cdfnames(j))) == 'edge_ref') then
            ncdf = ncdf - 1
          elseif (index(lcase(cdfnames(j)), 'nsd_') /= 0 .or. index(lcase(cdfnames(j)), 'nrc_') /= 0 .or. &
                  index(lcase(cdfnames(j)), 'nra_') /= 0) then
            ncdf = ncdf - 1
          elseif (allocated(fieldnames)) then !if -in option is present, discard unselected fields
            ffound = .false.
            do k = 1, size(fieldnames,1)
              if(trim(adjustl(cdfnames(j)))== trim(adjustl(fieldnames(k)))) ffound = .true.
            end do
            if (.not. ffound) ncdf = ncdf - 1
          end if
        end do
        if (allocated(cdfval)) deallocate(cdfval)
        allocate(cdfval(ncdf))
        ncdf = 0

        if(allocated(fieldnames)) then
          do k = 1, size(fieldnames,1)
            do j = 1, size(cdfnames,1)
              if(trim(adjustl(cdfnames(j))) == trim(adjustl(fieldnames(k))))  then
                ! Discard references
                if (trim(lcase(cdfnames(j))) == 'element_ref' .or. trim(lcase(cdfnames(j))) == 'face_ref' .or. &
                    trim(lcase(cdfnames(j))) == 'edge_ref') cycle
                if (index(lcase(cdfnames(j)), 'nsd_') /= 0 .or. index(lcase(cdfnames(j)),'nrc_') /= 0 .or. &
                    index(lcase(cdfnames(j)), 'nra_') /= 0) cycle
                if (vtk_var_xml_read('Cell', ncref, ncomp, trim(cdfnames(j)), temp, i) == 0) then
                  call info('    Reading cell field: '//trim(cdfnames(j)))
                  if (ncref == nc) then
                    ncdf = ncdf + 1
                    allocate(cdfval(ncdf)%val(ncomp,nc))
                    cdfval(ncdf)%val(:,:) = reshape(temp, [ncomp,nc])
                  else
                    call error('Wrong number of element values')
                  endif
                endif
                if(allocated(temp)) deallocate(temp)
              endif
            enddo
          enddo
        else
          do j = 1, size(cdfnames,1)
            ! Discard references
            if (trim(lcase(cdfnames(j))) == 'element_ref' .or. trim(lcase(cdfnames(j))) == 'face_ref' .or. &
                trim(lcase(cdfnames(j))) == 'edge_ref') cycle
            if (index(lcase(cdfnames(j)),'nsd_') /= 0 .or. index(lcase(cdfnames(j)),'nrc_') /= 0 .or. &
                index(lcase(cdfnames(j)), 'nra_') /= 0) cycle
            if (vtk_var_xml_read('Cell', ncref, ncomp, trim(cdfnames(j)), temp, i) == 0) then
              call info('    Reading cell field: '//trim(cdfnames(j)))
              if (ncref == nc) then
                ncdf = ncdf + 1
                allocate(cdfval(ncdf)%val(ncomp,nc))
                cdfval(ncdf)%val(:,:) = reshape(temp, [ncomp,nc])
              else
                call error('Wrong number of element values')
              endif
            endif
            if(allocated(temp)) deallocate(temp)
          enddo
        endif
      else
        ncdf = 0
      endif
    endif

    ! deallocate cdfval if no fields found
    if(allocated(cdfval) .and. ncdf == 0) then
      do j=1,size(cdfval,1)
        if(allocated(cdfval(j)%val)) deallocate(cdfval(j)%val)
      enddo
      deallocate(cdfval)
    endif

    ! Initialize celldata field names and number of parameters
    call suniqueI1P(cell_type, uct)
    ! Memory allocation: Element groups
    if(.not. allocated(pmh%pc(i)%el)) then
       allocate(pmh%pc(i)%el(size(uct,1)))
    elseif(size(pmh%pc(i)%el,1)<size(uct,1)) then
       call move_alloc(from=pmh%pc(i)%el, to=auxpmh%pc(1)%el)
       if(allocated(pmh%pc(i)%el)) deallocate(pmh%pc(i)%el)
       allocate(pmh%pc(i)%el(size(uct,1)))
       pmh%pc(i)%el(1:size(auxpmh%pc(1)%el,1)) = auxpmh%pc(1)%el(:)
       if(allocated(auxpmh%pc(1)%el)) deallocate(auxpmh%pc(1)%el)
    endif

    if(.not. allocated(auxpmh%pc(1)%el)) allocate(auxpmh%pc(1)%el(1))

    if(allocated(cdfval)) then
      do j=1, size(uct,1)
        ncdf = 0
        ! Memory allocation: Cell fields
        if(.not. allocated(pmh%pc(i)%el(j)%fi)) then
          allocate(pmh%pc(i)%el(j)%fi(size(cdfval,1)))
        elseif(size(pmh%pc(i)%el(j)%fi,1)<size(cdfval,1)) then
          call move_alloc(from=pmh%pc(i)%el(j)%fi, to=auxpmh%pc(1)%el(1)%fi)
          if(allocated(pmh%pc(i)%el(j)%fi)) deallocate(pmh%pc(i)%el(j)%fi)
          allocate(pmh%pc(i)%el(j)%fi(size(cdfval,1)))
          pmh%pc(i)%el(j)%fi(1:size(auxpmh%pc(1)%el(1)%fi,1)) = auxpmh%pc(1)%el(1)%fi(:)
          if(allocated(auxpmh%pc(1)%el(1)%fi)) deallocate(auxpmh%pc(1)%el(1)%fi)
        endif
        if(allocated(fieldnames)) then
          do k=1, size(fieldnames,1)
            do l=1,size(cdfnames,1)
              if (trim(lcase(cdfnames(j))) == 'element_ref' .or. trim(lcase(cdfnames(j))) == 'face_ref' .or. &
                  trim(lcase(cdfnames(j))) == 'edge_ref') cycle
              if (index(lcase(cdfnames(j)), 'nsd_') /= 0 .or. index(lcase(cdfnames(j)),'nrc_') /= 0 .or. &
                  index(lcase(cdfnames(j)), 'nra_') /= 0) cycle
              if(trim(adjustl(cdfnames(l))) == trim(adjustl(fieldnames(k)))) then
                ! if -in option is present, discard unselected fields
                ncdf = ncdf + 1
                pmh%pc(i)%el(j)%fi(ncdf)%name = trim(cdfnames(l))
                call set(pmh%pc(i)%el(j)%fi(ncdf)%param, p, np, fit=.true.)
              endif
            enddo
          enddo
        else
          do l = 1, size(cdfnames,1)
            if (trim(lcase(cdfnames(j))) == 'element_ref' .or. trim(lcase(cdfnames(j))) == 'face_ref' .or. &
                trim(lcase(cdfnames(j))) == 'edge_ref') cycle
            if (index(lcase(cdfnames(j)), 'nsd_') /= 0 .or. index(lcase(cdfnames(j)),'nrc_') /= 0 .or. &
                index(lcase(cdfnames(j)), 'nra_') /= 0) cycle
            ncdf = ncdf + 1
            pmh%pc(i)%el(j)%fi(ncdf)%name = trim(cdfnames(l))
            call set(pmh%pc(i)%el(j)%fi(ncdf)%param, p, np, fit=.true.)
          enddo
        endif
      enddo
      if(allocated(cdfnames)) deallocate(cdfnames)
    endif

    if(.not. allocated(auxpmh%pc(1)%el(1)%fi)) allocate(auxpmh%pc(1)%el(1)%fi(1))
    do j=1, size(uct,1)
      pmh%pc(i)%el(j)%type = vtk2pmhcelltype(int(uct(j)))
      lnn = FEDB(pmh%pc(i)%el(j)%type)%lnn
      call info('    Building group '//trim(string(j))//': '//trim(FEDB(pmh%pc(i)%el(j)%type)%desc))
      if(pmh%pc(i)%el(j)%type == 0) cycle
      nel = 0
      do k=1, nc
        ! Stores connectivities and references in PMH structure
        if(uct(j) == cell_type(k)) then
          nel = nel + 1
          call set(2, pmh%pc(i)%el(j)%nn, connect(offset(k)-lnn+1:offset(k))+1, nel, fit=tf)
          if(pmh%pc(i)%el(j)%type /= check_fe(.true.,1,1,0,0)) then
            call set(pmh%pc(i)%el(j)%ref, int(c_refs(k)), nel, fit=.true.)
          elseif(allocated(v_ref)) then
            call set(pmh%pc(i)%el(j)%ref, v_ref(connect(offset(k))+1), nel, fit=.true.)
          endif
        endif
        ! Saves VTU fields in a PMH structure. Manage memory allocation
        if(allocated(cdfval)) then
          do l=1,size(cdfval,1)
            if(.not. allocated(cdfval(l)%val)) cycle
            if(.not. allocated(pmh%pc(i)%el(j)%fi(l)%val)) then
              allocate(pmh%pc(i)%el(j)%fi(l)%val(size(cdfval(l)%val,1),size(cdfval(l)%val,2),1))
            elseif(size(pmh%pc(i)%el(j)%fi(l)%val,2) < nel) then
              if(allocated(auxpmh%pc(1)%el(1)%fi(1)%val)) deallocate(auxpmh%pc(1)%el(1)%fi(1)%val)
              call move_alloc(from=pmh%pc(i)%el(j)%fi(l)%val, to=auxpmh%pc(1)%el(1)%fi(1)%val)
              allocate(pmh%pc(i)%el(j)%fi(l)%val(size(cdfval(l)%val,1), nel,np))
              pmh%pc(i)%el(j)%fi(l)%val(1:size(auxpmh%pc(1)%el(1)%fi(1)%val,1), &
                                      & 1:size(auxpmh%pc(1)%el(1)%fi(1)%val,2), &
                                      & 1:size(auxpmh%pc(1)%el(1)%fi(1)%val,3)) = &
                                      & auxpmh%pc(1)%el(1)%fi(1)%val(:,:,:)
            elseif(size(pmh%pc(i)%el(j)%fi(l)%val,3) < np) then
              if(allocated(auxpmh%pc(1)%el(1)%fi(1)%val)) deallocate(auxpmh%pc(1)%el(1)%fi(1)%val)
              call move_alloc(from=pmh%pc(i)%el(j)%fi(l)%val, to=auxpmh%pc(1)%el(1)%fi(1)%val)
              allocate(pmh%pc(i)%el(j)%fi(l)%val(size(cdfval(l)%val,1), size(cdfval(l)%val,2),np))
              pmh%pc(i)%el(j)%fi(l)%val(1:size(auxpmh%pc(1)%el(1)%fi(1)%val,1), &
                                      & 1:size(auxpmh%pc(1)%el(1)%fi(1)%val,2), &
                                      & 1:size(auxpmh%pc(1)%el(1)%fi(1)%val,3)) = &
                                      & auxpmh%pc(1)%el(1)%fi(1)%val(:,:,:)
            endif
            pmh%pc(i)%el(j)%fi(l)%val(1:size(cdfval(l)%val,1),nel,np) = cdfval(l)%val(:,k)
          enddo
        endif
      enddo
      pmh%pc(i)%el(j)%nel = nel
      ! Reduce arrays to the right size
      call reduce(pmh%pc(i)%el(j)%nn,FEDB(pmh%pc(i)%el(j)%type)%lnn,pmh%pc(i)%el(j)%nel)
      call reduce(pmh%pc(i)%el(j)%ref,pmh%pc(i)%el(j)%nel)
      if(allocated(pmh%pc(i)%el(j)%fi)) then

        do l=1,size(pmh%pc(i)%el(j)%fi,1)
          if(allocated(temp2)) deallocate(temp2)
          allocate(temp2( &
                    & size(pmh%pc(i)%el(j)%fi(l)%val,1), &
                    & pmh%pc(i)%el(j)%nel , &
                    & size(pmh%pc(i)%el(j)%fi(l)%val,3) &
                 & ))
           temp2(:,1:pmh%pc(i)%el(j)%nel,1) = pmh%pc(i)%el(j)%fi(l)%val(:,1:pmh%pc(i)%el(j)%nel, np)
           call move_alloc(from=temp2, to=pmh%pc(i)%el(j)%fi(l)%val)
           if(allocated(temp2)) deallocate(temp2)
         enddo
       endif
    enddo

    ! Memory deallocation
    if(allocated(cdfval)) then
      do j=1,size(cdfval,1)
        if(allocated(cdfval(j)%val))deallocate(cdfval(j)%val)
      enddo
    endif

    if(allocated(X)) deallocate(X)
    if(allocated(Y)) deallocate(Y)
    if(allocated(Z)) deallocate(Z)
    if(allocated(connect)) deallocate(connect)
    if(allocated(offset)) deallocate(offset)
    if(allocated(cell_type)) deallocate(cell_type)
    if(allocated(uct)) deallocate(uct)
    if(allocated(v_ref)) deallocate(v_ref)
    if(allocated(c_refs)) deallocate(c_refs)
    if(allocated(aux_ref)) deallocate(aux_ref)
    if(allocated(pdfnames)) deallocate(pdfnames)
    if(allocated(cdfnames)) deallocate(cdfnames)
    if(allocated(pdfval)) deallocate(pdfval)
    if(allocated(cdfval)) deallocate(cdfval)
    if(allocated(temp)) deallocate(temp)
    if(allocated(temp2)) deallocate(temp2)

  enddo

  if (vtk_end_xml_read()/=0) stop 'Error'


end subroutine

!-----------------------------------------------------------------------
! suniqueI1P: returns the same values as in v with no repetitions.
! REMARK: res will not be sorted
!-----------------------------------------------------------------------
subroutine suniqueI1P(v, res)
integer(I1P), intent(in) :: v(:)
integer(I1P), allocatable, intent(out) :: res(:)
integer(I1P), dimension(size(v,1)) :: tmp
integer :: ios, n, i
character(maxpath) :: err_msg

n = 0
do i = 1, size(v,1) !loop to find the elements
  if (all(tmp(1:n) /= v(i))) then
    n = n + 1
    tmp(n) = v(i)
  endif
enddo
allocate(res(n), stat=ios, errmsg=err_msg)
if (ios /= 0) call error('unique, unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end subroutine

function vtk2pmhcelltype(vtktype) result(res)
  integer, intent(in) :: vtktype
  integer             :: res

  res = 0
  if(vtktype == 1) then ! VTK_VERTEX
    res = check_fe(.true., 1, 1, 0, 0)
  elseif(vtktype == 2) then ! VTK_POLY_VERTEX
    res = check_fe(.true., 1, 1, 0, 0)
  elseif(vtktype == 3) then ! VTK_LINE
    res = check_fe(.true., 2, 2, 1, 0)
  elseif(vtktype == 4) then ! VTK_POLY_LINE
    call error('Element type 4 (VTK_POLY_LINE) not supported!')
  elseif(vtktype == 5) then ! VTK_TRIANGLE
    res = check_fe(.true., 3, 3, 3, 0)
  elseif(vtktype == 6) then ! VTK_TRIANGLE_STRIP
    call error('Element type 6 (VTK_TRIANGLE_STRIP) not supported!')
  elseif(vtktype == 7) then ! VTK_POLYGON
    call error('Element type 7 (VTK_POLYGON) not supported!')
  elseif(vtktype == 8) then ! VTK_PIXEL
    res = check_fe(.true., 4, 4, 4, 0)
  elseif(vtktype == 9) then ! VTK_QUAD
    res = check_fe(.true., 4, 4, 4, 0)
  elseif(vtktype == 10) then ! VTK_TETRA
    res = check_fe(.true., 4, 4, 6, 4)
  elseif(vtktype == 11) then ! VTK_VOXEL
    res = check_fe(.true., 8, 8, 12, 6)
  elseif(vtktype == 12) then ! VTK_HEXAHEDRON
    res = check_fe(.true., 8, 8, 12, 6)
  elseif(vtktype == 13) then ! VTK_WEDGE
    res = check_fe(.true., 6, 6, 9, 5)
  elseif(vtktype == 14) then ! VTK_PYRAMID
    call error('Element type 14 (VTK_PYRAMID) not supported!')
  elseif(vtktype == 21) then ! VTK_QUADRATIC_EDGE
    res = check_fe(.false., 3, 2, 1, 0)
  elseif(vtktype == 22) then ! VTK_QUADRATIC_TRIANGLE
    res = check_fe(.false., 6, 3, 3, 0)
  elseif(vtktype == 23) then ! VTK_QUADRATIC_QUAD
    res = check_fe(.false., 8, 4, 4, 0)
  elseif(vtktype == 24) then ! VTK_QUADRATIC_TETRA
    res = check_fe(.false., 10, 4, 6, 4)
  elseif(vtktype == 25) then ! VTK_QUADRATIC_HEXAHEDRON
    res = check_fe(.false., 20, 8, 12, 6)
  endif

end function

!-----------------------------------------------------------------------
! type_cell: give the associated name of FE
!-----------------------------------------------------------------------
function pmh2vtkcelltype(pmhtype) result(res)
  integer, intent(in) :: pmhtype  !PMH cell type
  integer             :: res

  res = 0
  if(pmhtype == check_fe(.true., 1, 1, 0, 0)) then ! VTK_VERTEX
    res = 1
  elseif(pmhtype == check_fe(.true., 2, 2, 1, 0)) then ! VTK_LINE
    res = 3
  elseif(pmhtype == check_fe(.true., 3, 3, 3, 0)) then ! VTK_TRIANGLE
    res = 5
  elseif(pmhtype == check_fe(.true., 4, 4, 4, 0)) then ! VTK_QUAD
    res = 9
  elseif(pmhtype == check_fe(.true., 4, 4, 6, 4)) then ! VTK_TETRA
    res = 10
  elseif(pmhtype == check_fe(.true., 8, 8, 12, 6)) then ! VTK_HEXAHEDRON
    res = 12
  elseif(pmhtype == check_fe(.true., 6, 6, 9, 5)) then ! VTK_WEDGE
    res = 13
  elseif(pmhtype == check_fe(.false., 3, 2, 1, 0)) then ! VTK_QUADRATIC_EDGE
    res = 21
  elseif(pmhtype == check_fe(.false., 6, 3, 3, 0)) then ! VTK_QUADRATIC_TRIANGLE
    res = 22
  elseif(pmhtype == check_fe(.false., 8, 4, 4, 0)) then ! VTK_QUADRATIC_QUAD
    res = 23
  elseif(pmhtype == check_fe(.false., 10, 4, 6, 4)) then ! VTK_QUADRATIC_TETRA
    res = 24
  elseif(pmhtype == check_fe(.false., 20, 8, 12, 6)) then ! VTK_QUADRATIC_HEXAHEDRON
    res = 25
  else
    call error('Element type '//trim(string(pmhtype))//' not supported')
  endif

end function

!!-----------------------------------------------------------------------
!! PRIVATE PROCEDURES
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!! search_multiple: search the smallest value of 2 to the power of a that is bigger than b
!! 2^n*a > b  <=>  n > log2(b/a)
!!-----------------------------------------------------------------------
!integer function search_multiple(a,b)
!integer, intent(in) :: a, b
!
!if (b > a) then
!  search_multiple = int(2**real(ceiling(log(real(b)/a)/log(2.)))*a)
!else
!  search_multiple = a
!end if
!end function


end module
