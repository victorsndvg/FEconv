module module_pmh
!-----------------------------------------------------------------------
! Module to manage piecewise meshes (PMH)
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 19/01/2014
!
! PUBLIC PROCEDURES:
!   pmh2mfm: convert a PMH mesh into a MFM one
!
! REMARKS:
!   A mesh is divided into pieces 
!   Each piece has common vertices/nodes and it can contain one or several groups of elements
!   Each group of elements belong to a type of element defined in module_eltype
!   If either znod or z is known, the other variable can be constructed using el()%type, el()%nn, , el()%mm
!   For P1 elements, el()%nn and el()%mm are redundant
!
! PASOS A DAR:
!   Definición de pmh
!   Creacion de funciones para procesar las coordenadas, el nn, etc
!   Opcionalmente, creación de funciones para procesar refs, orientación, etc.
!   Conversor a/de mfm
!   Conversor a vtu
!   Uso en comsol
!
! MAS ACLARACIONES:
!   EL formato PMH sigue la misma regla que MFM para las conectividades, es decir:
!   1) mm siempre se guarda; 2) nn solo se guarda si no es P1; 3) z se define sobre los vértices
!   Por tanto, es responsabilidad del conversor a PMH el generar esta información (se creara
!   un módulo para que pueda ser llamado por todos los conversores, basado en las funciones de UNV)
!   las submallas no tienen por que estar ordenadas; es responsibilidad del lector de PMH de ordenarlas
!   para generar las matrices de referencias de MFM 
!-----------------------------------------------------------------------
use module_alloc
use module_args, only: get_post_arg
use module_fe_database_pmh
implicit none

!Constants
!EDGE_TRIA(i,j), vertex #i of edge #j of a triangle
integer, parameter :: EDGE_TRIA(2,3) = reshape([1,2, 2,3, 3,1], [2,3])

!EDGE_QUAD(i,j), vertex #i of edge #j of a quadrangle
integer, parameter :: EDGE_QUAD(2,4) = reshape([1,2, 2,3, 3,4, 4,1], [2,4])

!EDGE_TETRA(i,j), vertex #i of edge #j of a tetrahedron
integer, parameter :: EDGE_TETRA(2,6) = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4], [2,6])

!FACE_TETRA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: FACE_TETRA(3,4) = reshape([1,3,2, 1,4,3, 1,2,4, 2,3,4], [3,4])

!EDGE_HEXA(i,j), vertex #i of edge #j of a hexahedron
integer, parameter :: EDGE_HEXA(2,12) = reshape([1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5], [2,12])

!FACE_HEXA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: FACE_HEXA(4,6) = reshape([1,4,3,2, 1,5,8,4, 1,2,6,5, 5,6,7,8, 2,3,7,6, 3,4,8,7], [4,6])

!EDGE_WEDGE(i,j), vertex #i of edge #j of a wedge
integer, parameter :: EDGE_WEDGE(2,9) = reshape([1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4], [2,9])

!FACE_WEDGE(i,j), vertex #i of face #j of a wedge
integer, parameter :: FACE_WEDGE(4,5) = reshape([1,3,2,0, 1,4,6,3, 1,2,5,4, 4,5,6,0, 2,3,6,5], [4,5]) !note that faces 1 and 4 have only three valid vertices

!Types
type elgroup
  integer              :: type = 0 !element type (one of those defined in module_eltype)
  integer              :: nel  = 0 !total number of elements
  integer              :: lnn  = 0 !local number of nodes per element
  integer              :: lnv  = 0 !local number of vertices per element
  integer, allocatable :: nn(:,:)  !global numbering of nodes
  integer, allocatable :: mm(:,:)  !global numbering of vertices
  integer, allocatable :: ref(:)   !reference numbering
end type

type piece
  integer                    :: dim  = 0  ! space dimension of the node/vertex coordinates
  integer                    :: nnod = 0  ! total number of nodes
  integer                    :: nver = 0  ! total number of vertices
  real(real64),  allocatable :: znod(:,:) ! node coordinates
  real(real64),  allocatable :: z(:,:)    ! vertex coordinates
  type(elgroup), allocatable :: el(:)     ! element groups  
end type

type pmh_mesh
  type(piece), allocatable :: pc(:) ! pieces that compose the mesh
end type  

contains

!-----------------------------------------------------------------------
! pmh2mfm: convert a PMH mesh into a MFM one
!
! pmh is deallocated while MFM variables are being allocated
!-----------------------------------------------------------------------
subroutine pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
type(pmh_mesh), intent(inout) :: pmh
integer, intent(inout) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf 
integer, allocatable   :: nn(:,:), mm(:,:), nrv(:,:), nra(:,:), nrc(:,:), nsd(:)
real(real64), allocatable :: z(:,:)

integer :: i, ip, ig, pos, k, prev_nel, n, j, type_by_tdim(0:3), tmp_2d(2), tmp_3d(3)
integer, allocatable :: piece2save(:), ref(:,:), tmp_vf(:)
logical :: there_are_P1, there_are_other
character(maxpath) :: str

!check piece(s) to be saved
if (is_arg('-p')) then 
  str = get_post_arg('-p')
  if (str == '*') !save all pieces
    call alloc(piece2save, size(pmh%pc,1))
    piece2save = [(i, i=1, size(pmh%pc,1))]
  else !save a single piece, indicated in -p
    call alloc(piece2save, 1)
    piece2save(1) = int(str)
  end if
else !save only the first piece
  call alloc(piece2save, 1)
  piece2save(1) = 1
end if

!save a single piece !luego habrá que no diferenciar, ponerlo en un bucle...
if (size(piece2save,1) == 1) then
  ip = piece2save(1)

  !check whether all elements are P1 or not
  there_are_P1    = .false. !there are elements of type P1?
  there_are_other = .false. !there are elements of type different from P1?
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (FEDB(tp)%tdim > 0) then 
      there_are_P1    = there_are_P1    .or.      FEDB(tp)%nver_eq_nnod
      there_are_other = there_are_other .or. .not.FEDB(tp)%nver_eq_nnod
    end associate
  end do
  if (there_are_P1 .and. there_are_other) call error('(module_pmh/pmh2mfm) unable to convert P1 and non-P1 elements to MFM')
  if (there_are_other .and. (.not.allocated(pmh%pc(ip)%mm) .or. .not.allocated(pmh%pc(ip)%z))) call error('(module_pmh/pmh2mfm)'//&
  &' vertex coordinates are undefined for non-P1 elements, unable to convert to MFM')

  !check whether there is only one type of element for each topological dimension
  type_by_tdim = 0      !store the type of element for each topological dimension
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (tdim( FEDB(tp)%tdim ) == 0) then  
        type_by_tdim( FEDB(tp)%tdim ) = tp
      elseif (type_by_tdim( FEDB(tp)%tdim ) /= tp) then
        call error('(module_pmh/pmh2mfm) more that one type of element is defined for the same topological dimension: '//&
        &string(type_by_tdim(FEDB(tp)%tdim)//', '//string(tp)//'; unable to convert to MFM')
      end if
    end associate
  end do
  !maximal topological dimension
  max_tdim = 0 
  do i = 1, 3
    if (type_by_tdim(i) > 0) max_tdim  = i
  end do

  !save into MFM format
  nnod = pmh%pc(ip)%nnod
  nver = pmh%pc(ip)%nver
  lnn  = FEDB(type_by_tdim(max_tdim))%lnn
  lnv  = FEDB(type_by_tdim(max_tdim))%lnf
  lne  = FEDB(type_by_tdim(max_tdim))%lne
  lnf  = FEDB(type_by_tdim(max_tdim))%lnv

  !nel: total number of elements of maximal topological dimension
  nel = 0
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (FEDB( elg%type )%tdim == max_tdim) then
        nel = nel + elg%nel
      end if
    end associate
  end do

  !mm: concatenate vertex numbering of maximal topological dimension
  call alloc(mm, lnv, nel)
  prev_nel = 0
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (FEDB( elg%type )%tdim == max_tdim) then
        mm(1:lnv, prev_nel+1:prev_nel+elg%nel) = elg%mm
        call dealloc(elg%mm)
        prev_nel = prev_nel + elg%nel
      end if
    end associate
  end do

  !nsd: concatenate references of maximal topological dimension
  call alloc(nsd, nel)
  prev_nel = 0
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (FEDB( elg%type )%tdim == max_tdim) then
        nsd(prev_nel+1:prev_nel+elg%nel) = ref
        call dealloc(elg%nsd)
        prev_nel = prev_nel + elg%nel
      end if
    end associate
  end do

  !nn: concatenate node numbering of maximal topological dimension
  if (there_are_other) then
    call alloc(nn, lnn, nel)
    prev_nel = 0
    do ig = 1, size(pmh%pc(ip)%el, 1)
      associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
        if (FEDB( elg%type )%tdim == max_tdim) then
          nn(1:lnn, prev_nel+1:prev_nel+elg%nel) = elg%nn
          call dealloc(elg%nn)
          prev_nel = prev_nel + elg%nel
        end if
      end associate
    end do
  end if

  !nrv: visit element groups of tdim = 0 to set vertex references
  if (max_tdim > 0) then
    !STEP 1: create ref to collect vertices and references stored in PMH, orderly
    n = 0
    do ig = 1, size(pmh%pc(ip)%el, 1)
      associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
        if (FEDB( elg%type )%tdim == 0) then
          do k = 1, elg%nel
            tmp_2d = [elg%mm(1,k), elg%ref(k)]
            call insert_row_sorted(ref, tmp_2d, used=n, fit=[.false.,.true.])
          end do
          call dealloc(elg%ref)
        end if
      end associate
    end do
    call reduce(ref, n, 2)
    !STEP2: for every vertex in mm, check whether it is in ref
    do k = 1, nel
      do j = 1, FB(tp_max_tdim)%lnv
        pos = find_row_sorted(ref(1,:), mm(j,k), n)
        if (pos > 0) nrv(j,k) = ref(pos, 2)
      end do
    end do
    call dealloc(ref)
  end if

  !nra: visit element groups of tdim = 1 to set edge references
  if (max_tdim > 1) then
    !STEP 1: create ref to collect edge vertices and references stored in PMH, orderly
    n = 0
    do ig = 1, size(pmh%pc(ip)%el, 1)
      associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
        if (FEDB( elg%type )%tdim == 1) then
          do k = 1, elg%nel
            tmp_3d = [sort(elg%mm(1:2,k)), elg%ref(k)]
            call insert_row_sorted(ref, tmp_3d, used=n, fit=[.false.,.true.])
          end do
          call dealloc(elg%ref)          
        end if
      end associate
    end do
    call reduce(ref, n, 3)
    !STEP2: for every edge in mm, check whether it is in ref
    do k = 1, nel
      do j = 1, FB(tp_max_tdim)%lne
        tmp_2d = sort(mm(FB(tp_max_tdim)%edge(:,j),k))
        pos = find_row_sorted(ref(1:2,:), tmp_2d, n)
        if (pos > 0) nra(j,k) = ref(pos, 3)
      end do
    end do
    call dealloc(ref)    
  end if

  !nrc: visit element groups of tdim = 2 to set edge references
  if (max_tdim > 2) then
    !STEP 1: create ref2 to collect face vertices and references stored in PMH, orderly
    n = 0
    associate(v_f => FB(tp_max_tdim)%lnv_f) !v_f: vertices per face in the max. tdim element
      call alloc(tmp_vf, v_f+1)
      do ig = 1, size(pmh%pc(ip)%el, 1)
        associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
          if (FEDB( elg%type )%tdim == 2 .and. FEDB(elg%type)%lnv == FB(tp_max_tdim)%lnv_f) then
            do k = 1, elg%nel
              tmp_vf = [sort(elg%mm(1:v_f, k)), elg%ref(k)]
              call insert_row_sorted(ref, tmp_vf, used=n, fit=[.false.,.true.])
            end do
            call dealloc(elg%ref)
          end if
        end associate
      end do
      call reduce(ref, n, v_f+1)
      !STEP2: for every face in mm, check whether it is in ref
      call alloc(tmp_vf, v_f)
      do k = 1, nel
        do j = 1, FB(tp_max_tdim)%lnf
          tmp_vf = sort(mm(FB(tp_max_tdim)%face(:,j),k))
          pos = find_row_sorted(ref(1:v_f,:), tmp_vf, n)
          if (pos > 0) nrc(j,k) = ref(pos, v_f+1)
        end do
      end do
    end associate  
    call dealloc(ref)
  end if

  !z: save vertex coordinates
  call alloc(z, FEDB(tp_max_tdim)%dim, pmh%pc(ip)%nver)
  z = pmh%pc(ip)%z
  call dealloc(pmh%pc(ip)%z)

else !there are several pieces to save
  !if (is_arg('-glue')) then consider to identify nodes/vertices in the boundary of two pieces
  call error('(module_pmh/pmh2mfm) saving several pieces is not implemented')
end if

!hay que hacer subrutinas para dado p2, construir P1 en todos los elementos)
end subroutine

end module