module module_vtu
!-----------------------------------------------------------------------
! Module to manage VTU files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 15/05/2013
!
! PUBLIC PROCEDURES:
!   save_vtu: saves a mesh in a VTU format file
!   type_cell: give the associated name of FE
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_report, only:error, report_option
use module_files, only:get_unit
use module_convers, only:string
use module_alloc
use module_set, only: sunique
use LIB_VTK_IO
use LIB_VTK_IO_READ
use module_writevtu
use module_pmh
use module_fe_database_pmh
implicit none

!Constants
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
private :: true64, save_w_field !save_w_field must change

contains

!-----------------------------------------------------------------------
! save_vtu: save mesh with references
!-----------------------------------------------------------------------
subroutine save_vtu(filename, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
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



subroutine save_vtu2(filename, pmh)
  character(len=*), intent(in)  :: filename
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  real(real64), allocatable     :: znod(:,:)
  integer, allocatable          :: connect(:), offset(:), celltypes(:)
  integer, allocatable          :: v_ref(:), e_ref(:), f_ref(:), el_ref(:)
  integer :: i, j, k, l, m
  integer :: nnod, nel, conlen, tp, lnn, lnv, tnvpc, maxtopdim
  logical :: all_P1

  call VTU_open(trim(filename)) 
  ! Loop in pieces
  do i=1,size(pmh%pc,1)
    nel = 0 ! Number of VTK elements
    tnvpc = 0 ! Total number of vertex per cell
    maxtopdim = 0
    call build_node_coordinates(pmh%pc(i), i, all_P1, znod)
    if(allocated(v_ref)) deallocate(v_ref)
    if(.not. all_P1)then
      allocate(v_ref(size(znod,2)))
    else
      allocate(v_ref(pmh%pc(i)%nver))
    endif
    v_ref = 0
    ! Calc max topological dimension
    do j=1,size(pmh%pc(i)%el,1)
      maxtopdim = max(maxtopdim,FEDB(pmh%pc(i)%el(j)%type)%tdim)
    enddo
    ! Loop in element groups
    do j=1,size(pmh%pc(i)%el,1)
      tp = pmh%pc(i)%el(j)%type
      print*, 'Element type: '//trim(FEDB(tp)%desc)
      if(.not. FEDB(tp)%nver_eq_nnod) then
          lnn = FEDB(tp)%lnn
          ! Build vtk conectivity array
          call set(connect, pack(reshape(pmh%pc(i)%el(j)%nn-1, (/1,pmh%pc(i)%el(j)%nel*lnn/)),.true.), &
            &(/(k,k=tnvpc+1,tnvpc+pmh%pc(i)%el(j)%nel*lnn)/), fit=.false.)
          ! Build vtk offset array
          call set(offset, (/(tnvpc+k*lnn, k=1,pmh%pc(i)%el(j)%nel)/) , &
            &(/(k,k=nel,nel+pmh%pc(i)%el(j)%nel)/), fit=.false.)
          tnvpc = tnvpc + pmh%pc(i)%el(j)%nel*lnn-1+lnv
      else
          lnv = FEDB(tp)%lnv
          ! Build vtk conectivity array
          call set(connect, pack(reshape(pmh%pc(i)%el(j)%mm-1, (/1,pmh%pc(i)%el(j)%nel*lnv/)),.true.), &
            &(/(k,k=tnvpc+1,tnvpc+pmh%pc(i)%el(j)%nel*lnv)/), fit=.false.)
          ! Build vtk offset array
          call set(offset, (/(tnvpc+k*lnv, k=1,pmh%pc(i)%el(j)%nel)/) , &
            &(/(k,k=nel+1,nel+1+pmh%pc(i)%el(j)%nel)/), fit=.false.)
          tnvpc = tnvpc + pmh%pc(i)%el(j)%nel*lnv
      endif
      ! Build vtk cell types array
      call set(celltypes, (/(pmh2vtkcelltype(tp), k=nel,nel+pmh%pc(i)%el(j)%nel)/) , &
        &(/(k,k=nel+1,nel+1+pmh%pc(i)%el(j)%nel)/), fit=.false.)
      if(FEDB(tp)%tdim == maxtopdim) then
        ! 'Building subdomain references array ... '
        call set(el_ref, (/(pmh%pc(i)%el(j)%ref(k), k=1,pmh%pc(i)%el(j)%nel)/) , &
          &(/(k,k=nel+1,nel+1+pmh%pc(i)%el(j)%nel)/), fit=.false.)
      elseif(FEDB(tp)%tdim == 2) then
        ! 'Building face references array ... '
        call set(f_ref, (/(pmh%pc(i)%el(j)%ref(k), k=1,pmh%pc(i)%el(j)%nel)/) , &
          &(/(k,k=nel+1,nel+1+pmh%pc(i)%el(j)%nel)/), fit=.false.)
      elseif(FEDB(tp)%tdim == 1) then
        ! 'Building edge references array ... '
        call set(e_ref, (/(pmh%pc(i)%el(j)%ref(k), k=1,pmh%pc(i)%el(j)%nel)/) , &
          &(/(k,k=nel+1,nel+1+pmh%pc(i)%el(j)%nel)/), fit=.false.)
      elseif(FEDB(tp)%tdim == 0) then
        ! 'Building vertex references array ... '
        do k=1,pmh%pc(i)%el(j)%nel
          v_ref(pmh%pc(i)%el(j)%mm(:,k)) = pmh%pc(i)%el(j)%ref(k)
        enddo
        ! Vertex as celldata
!        call set(v_ref, (/(pmh%pc(i)%el(j)%ref(k), k=1,pmh%pc(i)%el(j)%nel)/) , &
!          &(/(k,k=nel+1,nel+1+pmh%pc(i)%el(j)%nel)/), fit=.false.)
      endif
      nel = nel + pmh%pc(i)%el(j)%nel
    enddo
    if(allocated(connect)) call reduce(connect, tnvpc)
    if(allocated(offset)) call reduce(offset, nel)
    if(allocated(celltypes)) call reduce(celltypes, nel)
    if(pmh%pc(i)%dim < 3) then
      if(.not. all_P1) then
        nnod = size(znod,2)
        if (vtk_geo_xml(nnod, nel, pack(znod(1,:),.true.), pack(znod(2,:),.true.), (/(0._real64,l=1,nnod)/)) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      else
        nnod = pmh%pc(i)%nver
        if (vtk_geo_xml(pmh%pc(i)%nver, nel, &
        &  pack(pmh%pc(i)%z(1,:),.true.), pack(pmh%pc(i)%z(2,:),.true.), (/(0._real64,l=1,pmh%pc(i)%nver)/)) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      endif
    else
      if(.not. all_P1) then
        nnod = size(znod,2)
        if (vtk_geo_xml(nnod, nel, pack(znod(1,:),.true.), pack(znod(2,:),.true.), pack(znod(3,:),.true.)) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      else
        nnod = pmh%pc(i)%nver
        if (vtk_geo_xml(pmh%pc(i)%nver, nel, &
        &  pack(pmh%pc(i)%z(1,:),.true.), pack(pmh%pc(i)%z(2,:),.true.), pack(pmh%pc(i)%z(3,:),.true.)) /= 0) &
        &  stop 'writeVTU: vtk_geo_xml error 1'
      endif
    endif
    if(vtk_con_xml(nel,connect,offset,int(celltypes,I1P)) /= 0) stop 'writeVTU: vtk_con_xml error 1'

    call VTU_begin_pointdata()
    if(allocated(v_ref) .and. size(pack(v_ref,v_ref/=0))>0) then
      Print*, '  Vertex references found!'
      call reduce(v_ref, nnod)
      if(vtk_var_xml(nnod, 'vertex_ref', v_ref) /= 0) stop 
    endif
    call VTU_end_pointdata()
    call VTU_begin_celldata()
    if(allocated(el_ref) .and. size(pack(el_ref,el_ref/=0))>0) then
      Print*, '  Subdomain references found!'
      call reduce(el_ref, nel)
      if(vtk_var_xml(nel, 'element_ref', el_ref) /= 0) stop 
    endif
    if(allocated(f_ref)) then
      Print*, '  Face references found!'
      call reduce(f_ref, nel)
      if(vtk_var_xml(nel, 'face_ref', f_ref) /= 0) stop 
    endif
    if(allocated(e_ref)) then
      Print*, '  Edge references found!'
      call reduce(e_ref, nel)
      if(vtk_var_xml(nel, 'edge_ref', e_ref) /= 0) stop 
    endif
    call VTU_end_celldata()
    if(vtk_geo_xml() /= 0) stop 'writeVTU: vtk_geo_xml_closep error 1'
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

subroutine load_vtu(filename, pmh)
  character(len=*), intent(in)          :: filename
  type(pmh_mesh), intent(inout)         :: pmh

  ! Inital settings
  call report_option('level', 'stdout')
 
  ! Open VTU file and reads the mesh
  call info('Reading VTU file ...')
  call read_vtu(filename, pmh)

  call build_vertices(pmh)

end subroutine


!-----------------------------------------------------------------------
! read_VTU(filename, pmh): read VTU file
!-----------------------------------------------------------------------
! this:   vtu type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine read_vtu(filename, pmh)

  character(len=*), intent(in)  :: filename
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer(I4P)              :: np, nn, nc, nco, nnref, ncref, ncomp
  real(R8P), allocatable    :: X(:), Y(:), Z(:)
  integer, allocatable      :: v_ref(:), aux_ref(:), c_refs(:)
  integer(I4P), allocatable :: connect(:), offset(:)
  integer(I1P), allocatable :: cell_type(:),uct(:)
  integer                   :: i,j,k, lnn, ini, fin


  ! Reads the VTU file header and allocates the number of pieces

  print*, 'Reading VTK header ...'
  if (vtk_ini_xml_read('Binary',filename,'UnstructuredGrid', np)/=0) stop 'Error'
  print*, 'Pieces founded: '//trim(string(np))

  if(allocated(pmh%pc)) deallocate(pmh%pc)
  allocate(pmh%pc(np))

  do i=1, np
    print*, 'Piece '//trim(string(i))//':'
    print*, '  Reading coordinates ...'
    if (vtk_geo_xml_read(nn, nc, X, Y, Z, i)/=0) stop 'Error'
    pmh%pc(i)%nnod = nn
    pmh%pc(i)%dim = 3
    if(allocated(pmh%pc(i)%z)) deallocate(pmh%pc(i)%z)
    allocate(pmh%pc(i)%z(3,nn))

    
    call set_row(pmh%pc(i)%z, X, 1, fit=[.true.,.true.])
    call set_row(pmh%pc(i)%z, Y, 2, fit=[.true.,.true.])
    call set_row(pmh%pc(i)%z, Z, 3, fit=[.true.,.true.])


    print*, '  Reading conectivities ... '
    if (vtk_con_xml_read(nc,connect,offset,cell_type,i)/=0) stop 'Error'

    ! Total array of cell references (celldata)
    if(allocated(c_refs)) deallocate(c_refs)
    allocate(c_refs(nc))
    c_refs = 0

    print*, '  Reading references ... '
    ! VTK+ vertex references (pointdata)
    if(vtk_var_xml_read('Node',nnref,ncomp,'vertex_ref',v_ref,i) == 0 ) then
      if((nnref/=nn .or. ncomp /= 1) .and. allocated(v_ref)) deallocate(v_ref)
    endif
    ! VTK+ edge references (celldata)
    if(vtk_var_xml_read('Cell',ncref,ncomp,'edge_ref',aux_ref,i) == 0) then
      if((ncref/=nc .or. ncomp/= 1) .and. allocated(aux_ref)) then
        deallocate(aux_ref)
      elseif(allocated(aux_ref)) then
         call add(c_refs,aux_ref,(/(j,j=1,nc)/),fit=.true.)
         deallocate(aux_ref)
      endif
    endif
    ! VTK+ face references (celldata)
    if(vtk_var_xml_read('Cell',ncref,ncomp,'face_ref',aux_ref,i) == 0) then
      if((ncref/=nc .or. ncomp/=1) .and. allocated(aux_ref)) then
        deallocate(aux_ref)
      elseif(allocated(aux_ref)) then
        call add(c_refs,aux_ref,(/(j,j=1,nc)/),fit=.true.)
        deallocate(aux_ref)
      endif
    endif
    ! VTK+ element references (celldata)
    if(vtk_var_xml_read('Cell',ncref,ncomp,'element_ref',aux_ref,i) == 0) then
      if((ncref/=nc .or. ncomp /= 1) .and. allocated(aux_ref)) then 
        deallocate(aux_ref)
      elseif(allocated(aux_ref)) then
        call add(c_refs,aux_ref,(/(j,j=1,nc)/),fit=.true.)
        deallocate(aux_ref)
      endif
    endif


    uct = uniqueI1P(cell_type)
    if(allocated(pmh%pc(i)%el)) deallocate(pmh%pc(i)%el)
!    if(allocated(v_ref)) then 
!      allocate(pmh%pc(i)%el(size(uct,1)+1))
!    else
      allocate(pmh%pc(i)%el(size(uct,1)))
!    endif
    do j=1, size(uct,1)
      pmh%pc(i)%el(j)%type = vtk2pmhcelltype(int(uct(j)))
      lnn = FEDB(pmh%pc(i)%el(j)%type)%lnn
      print*, '  Building group '//trim(string(j))//': '//trim(FEDB(pmh%pc(i)%el(j)%type)%desc)
      if(pmh%pc(i)%el(j)%type == 0) cycle
      do k=1, nc
        if(uct(j) == cell_type(k)) then
          pmh%pc(i)%el(j)%nel = pmh%pc(i)%el(j)%nel + 1
          call set_col(pmh%pc(i)%el(j)%nn, connect(offset(k)-lnn+1:offset(k))+1, pmh%pc(i)%el(j)%nel, fit=[.true.,.false.])
          if(pmh%pc(i)%el(j)%type /= check_fe(.true.,1,1,0,0)) then
            call set(pmh%pc(i)%el(j)%ref, int(c_refs(k)), pmh%pc(i)%el(j)%nel, fit=.true.)
          elseif(allocated(v_ref)) then
            call set(pmh%pc(i)%el(j)%ref, v_ref(connect(offset(k))+1), pmh%pc(i)%el(j)%nel, fit=.true.)
          endif
        endif
      enddo
      call reduce(pmh%pc(i)%el(j)%nn,FEDB(pmh%pc(i)%el(j)%type)%lnn,pmh%pc(i)%el(j)%nel)
      call reduce(pmh%pc(i)%el(j)%ref,pmh%pc(i)%el(j)%nel)
    enddo

!    if(allocated(v_ref)) then 
!      j = size(uct,1) + 1
!      pmh%pc(i)%el(j)%type = check_fe(.true.,1,1,0,0)
!      call set_row(pmh%pc(i)%el(j)%nn, (/(j,j=1,nn)/), 1, fit=[.true.,.true.])
!print*, 'nn',trim(string(pack(pmh%pc(i)%el(j)%nn,.true.)))
!      call set(pmh%pc(i)%el(j)%ref, int(v_ref), (/(j,j=1,nn)/), fit=.true.)
!print*, 'ref',trim(string(pmh%pc(i)%el(j)%ref))
!    endif

  enddo

    if (vtk_end_xml_read()/=0) stop 'Error'


end subroutine

!-----------------------------------------------------------------------
! unique (function): returns the same values as in v with no repetitions.
! REMARK: res will not be sorted
!-----------------------------------------------------------------------
function uniqueI1P(v) result(res)
integer(I1P), intent(in) :: v(:)
integer(I1P), allocatable :: res(:)
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
end function

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



end module
