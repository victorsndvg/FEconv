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
use module_report, only:error
use module_convers, only:string
use module_alloc_int_r1, only: ssort
use module_set, only: sunique
use module_writevtu
implicit none

!Constants
!edge_tria(i,j), vertice i de la arista j de un triangulo
integer, parameter, dimension(2,3) :: edge_tria = reshape([1,2, 2,3, 3,1], [2,3])
!edge_tetra(i,j), vertice i de la arista j de un tetraedro
integer, parameter, dimension(2,6) :: edge_tetra = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4], [2,6])
!face_tetra(i,j), vertice i de la cara j de un tetraedro
integer, parameter, dimension(3,4) :: face_tetra = reshape([1,3,2, 1,4,3, 1,2,4, 2,3,4], [3,4])

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
      elv = 0._real64
      do k = 1, nel
        do l = 1, lne; do m = 1, 2
          elv(nn(edge_tetra(m,l),k)) = max(elv(nn(edge_tetra(m,l),k)), true64(nra(l,k)==ref(j)))
        end do; end do
        do l = 1, lnn-lnv
          eln(nn(l+lnv,k)) = max(eln(nn(l+lnv,k)), min(eln(nn(edge_tetra(1,l),k)),eln(nn(edge_tetra(2,l),k))))
        end do
      end do
      call VTU_write_pointdata(elv, 'nra_'//trim(string(ref(j))), 'scalar')
    enddo
  endif
  !nrc
  call sunique(pack(nrc, nrc/=0), ref)
  if (size(ref, 1) > 0) then 
    call ssort(ref)
    do j = 1, size(ref, 1)
      elv = 0._real64
      do k = 1,nel
        do l = 1, lnf; do m = 1,3
          elv(nn(face_tetra(m,l),k)) = max(elv(nn(face_tetra(m,l),k)), true64(nrc(l,k)==ref(j)))
        enddo; enddo
        do l = 1, lnn-lnv
          eln(nn(l+lnv,k)) = max(eln(nn(l+lnv,k)), min(eln(nn(face_tetra(1,l),k)),eln(nn(face_tetra(2,l),k)),&
          eln(nn(face_tetra(3,l),k))))
        end do
      enddo
      call VTU_write_pointdata(eln, 'nrc_'//trim(string(ref(j))), 'scalar')
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

end module
