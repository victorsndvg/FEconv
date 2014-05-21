module module_freefem
!-----------------------------------------------------------------------
! Module to manage MESH (FreeFem++) meshes (PMH)
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 20/05/2014
!
! PUBLIC PROCEDURES:
! save_freefem: save a PMH structure into a MESH (FreeFem++) file
!-----------------------------------------------------------------------
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_convers, only: string, int
use module_alloc, only: alloc, dealloc 
use module_args, only: is_arg, get_post_arg
use module_feed, only: feed, empty
use module_fe_database_pmh, only: FEDB, check_fe
use module_pmh, only: pmh_mesh
implicit none

contains

!-----------------------------------------------------------------------
! save_freefem: save a PMH structure into a MESH (FreeFem++) file
!
! pmh is deallocated while variables are being saved
!-----------------------------------------------------------------------
subroutine save_freefem(outfile, iu, pmh)
character(*),   intent(in)    :: outfile
integer,        intent(in)    :: iu
type(pmh_mesh), intent(inout) :: pmh

integer :: i, ipp, ip, ig, k, j, type_by_tdim(0:3), prev_max_tdim, res, max_tdim, ios, nel, ntri, nver, dim
integer, allocatable :: piece2save(:), tet_piece(:), tri_piece(:), nver_piece(:)
character(maxpath) :: str, cad

!check piece(s) to be saved
if (is_arg('-p')) then !save a single piece, indicated after -p
  str = get_post_arg('-p')
  call alloc(piece2save, 1)
  piece2save(1) = int(str)
else !save all pieces
  call alloc(piece2save, size(pmh%pc,1))
  piece2save = [(i, i=1, size(pmh%pc,1))]
end if
if (is_arg('-glue')) then 
  call info('(module_freefem/save_freefem) option -glue not implemented yet')
end if

!testing and calculation of max_tdim
type_by_tdim  = 0 !store the type of element for each topological dimension
prev_max_tdim = 0 !store the maximal topological dimension
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_freefem/save_freefem) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (tp /= check_fe(.true., 4, 4, 6, 4) .and. tp /= check_fe(.true., 3, 3, 3, 0)) then
        call info('(module_freefem/save_freefem) element type '//trim(FEDB(tp)%desc)//' found; those elements cannot be saved'//&
        &' in FreeFem++ format and they will be discarded')
        cycle
      end if  
      !check whether there is only one type of element for each topological dimension
      if (type_by_tdim( FEDB(tp)%tdim ) == 0) then  
        type_by_tdim( FEDB(tp)%tdim ) = tp
      elseif (type_by_tdim( FEDB(tp)%tdim ) /= tp) then
        call error('(module_freefem/save_freefem) more that one type of element is defined for the same topological dimension: '//&
        &string(type_by_tdim(FEDB(tp)%tdim))//', '//string(tp)//'; unable to convert to FreeFem++')
      end if
    end associate
  end do
  max_tdim = 0 
  do i = 1, 3
    if (type_by_tdim(i) > 0) max_tdim  = i
  end do
  if (prev_max_tdim == 0) then  
    prev_max_tdim = max_tdim
  elseif (prev_max_tdim /= max_tdim) then
    call error('(module_freefem/save_freefem) there are pieces with different maximal topological dimension; unable to convert '//&
    &'to FreeFem++')
  end if
end do
if (max_tdim < 3) call error('(module_freefem/save_freefem) only the conversion to a thetrahedral FreeFem++ mesh is implemented.')
!testing and calculation of dim
dim = -1
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  !check whether there is only one coordinates dimension for all pieces
  if (dim == -1) then  
    dim = pmh%pc(ip)%dim
  elseif (pmh%pc(ip)%dim /= dim) then
    call error('(module_freefem/save_freefem) different coordinates dimensions in different pieces: '//&
    &string(pmh%pc(ip)%dim)//', '//string(dim)//'; unable to convert to MFM')
  end if
end do

!store variables nnod, nver, nel for selected pieces
if (allocated(tet_piece)) deallocate(tet_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(tet_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable tet_piece: '//trim(cad))
if (allocated(tri_piece)) deallocate(tri_piece); !nel_piece(ipp):  global numbering for the last element of piece #ipp
allocate(tri_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable tri_piece: '//trim(cad))
if (allocated(nver_piece)) deallocate(nver_piece); !nver_piece(ipp): global numbering for the last vertex of piece #ipp
allocate(nver_piece(0:size(piece2save,1)), stat = res, errmsg = cad)
if (res /= 0) call error('(module_freefem/save_freefem) Unable to allocate variable nver_piece: '//trim(cad))
tet_piece(0) = 0; tri_piece(0) = 0; nver_piece(0) = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  nver_piece(ipp) = nver_piece(ipp-1) + pmh%pc(ipp)%nver
  tet_piece(ipp)  =  tet_piece(ipp-1) 
  tri_piece(ipp)  =  tri_piece(ipp-1) 
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == check_fe(.true., 4, 4, 6, 4)) then
         tet_piece(ipp) =  tet_piece(ipp) + elg%nel
      elseif (elg%type == check_fe(.true., 3, 3, 3, 0)) then
         tri_piece(ipp) =  tri_piece(ipp) + elg%nel
      end if
    end associate
  end do
end do
nel  =  tet_piece(size(piece2save,1))
ntri =  tri_piece(size(piece2save,1))
nver = nver_piece(size(piece2save,1))

!store variables nver, nel, ntri for selected pieces
open (unit=iu, file=outfile, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('save/open, #'//trim(string(ios)))
call feed(iu, string(nver)); call feed(iu, string(nel)); call feed(iu, string(ntri));  call empty(iu)

!write(iu, '(a/)') 'MeshVersionFormatted 1')
!write(iu, '(a)')  'Vertices'
!write(iu, '(a)')  trim(string(nver))

!store vertex coordinates and reference 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do j = 1, pmh%pc(ip)%nver
    do i = 1, pmh%pc(ip)%dim; call feed(iu, string(pmh%pc(ip)%z(i,j))); end do
    do i = 1, 3-dim;          call feed(iu, string(0.));                end do
    call feed(iu, string(0))
    call empty(iu)
  end do
end do

!store tet conectivities and references
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == check_fe(.true., 4, 4, 6, 4)) then
        do k = 1, elg%nel
          do i = 1, 4; call feed(iu, string(nver_piece(ipp-1)+elg%mm(i,k))); end do
          call feed(iu, string(elg%ref(k)))
          call empty(iu)
        end do
      end if
    end associate
  end do
end do

!store tri conectivities and references
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (elg%type == check_fe(.true., 3, 3, 3, 0)) then
        do k = 1, elg%nel
          do i = 1, 3; call feed(iu, string(nver_piece(ipp-1)+elg%mm(i,k))); end do
          call feed(iu, string(elg%ref(k)))
          call empty(iu)
        end do
      end if
    end associate
  end do
end do

close(iu)
end subroutine

end module
