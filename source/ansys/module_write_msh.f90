module module_write_msh_fcnv

!-----------------------------------------------------------------------
! Module to manage MSH (Ansys) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
!-----------------------------------------------------------------------
use basicmod
use module_transform_fcnv, only: to_l1
use module_pmh_fcnv
use module_utils_msh_fcnv
implicit none

contains

subroutine write_msh_header(iu, pmh, maxtopdim)
  integer, intent(in)           :: iu ! file descriptor
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer,        intent(inout) :: maxtopdim
  integer                       :: ip, ig, ios, totnver, totncells, totnfaces
  logical                       :: all_P1
  character(len=MAXPATH)        :: auxch1,auxch2,auxheader

  maxtopdim = 0
  totnver = 0
  totncells = 0
  totnfaces = 0
  all_P1 = .false.

  do ip=1,size(pmh%pc)
    do ig = 1, size(pmh%pc(ip)%el,1)
      if (FEDB(pmh%pc(ip)%el(ig)%type)%tdim > maxtopdim) then
        maxtopdim = FEDB(pmh%pc(ip)%el(ig)%type)%tdim
        totncells = pmh%pc(ip)%el(ig)%nel
      elseif (FEDB(pmh%pc(ip)%el(ig)%type)%tdim == maxtopdim) then
        totncells = totncells + pmh%pc(ip)%el(ig)%nel
      end if

      if (.not. FEDB(pmh%pc(ip)%el(ig)%type)%nver_eq_nnod) then
        call error('Must be a Lagrange P1 mesh...')
      end if

    end do

  totnver = totnver + pmh%pc(ip)%nver

  end do

  print*, 'Writing MSH headers'

  write(unit=iu, fmt='(a)', iostat = ios)  '(0 grid written by FEconv '
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(unit=iu, fmt='(a)', iostat = ios)  '   nodes:       (10 (id start end type) (x y z ...))'
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(unit=iu, fmt='(a)', iostat = ios)  '   faces:       (13 (id start end type elemType)'
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(unit=iu, fmt='(a)', iostat = ios)  '                (v-0 v-1 .. v-n right-cell left-cell ...))'
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(unit=iu, fmt='(a)', iostat = ios)  '   cells:       (12 (id start end type elemtype))'
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(unit=iu, fmt='(a)', iostat = ios)  '   zones:       (45 (id type name)())'
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(unit=iu, fmt='(a)', iostat = ios)  ')'
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))


  write(auxch1,'(I1)') maxtopdim
  write(auxheader,*) '(2 ',trim(adjustl(auxch1)),')'
  write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(auxheader))
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(auxch1,'(z20)') 1
  write(auxch2,'(z20)') totnver
  write(auxheader,*) '(10 (0 ',trim(adjustl(auxch1)),' ', trim(adjustl(auxch2)), ' 0))'
  write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(auxheader))
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))
  write(auxch2,'(z20)') totncells
  write(auxheader,*) '(12 (0 ',trim(adjustl(auxch1)),' ', trim(adjustl(auxch2)), ' 0))'
  write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(auxheader))
  if (ios /= 0) call error('module_write_msh/write_header # '//trim(string(ios)))

end subroutine


subroutine write_msh_nodes(iu, pmh, maxref,z)
  integer,               intent(in) :: iu ! file descriptor
  type(pmh_mesh),     intent(inout) :: pmh ! pmh_mesh
  integer,            intent(inout) :: maxref
  real(real64), allocatable, optional, intent(out) :: z(:,:)
  integer                       :: ip, i, ios, totnver, maxdim, ini, cnt
  character(len=MAXPATH)        :: ini_ch,end_ch,zone_ch,dim_ch, aux
  integer, allocatable          :: totrefs(:)


  totnver = 0
  maxdim = 0
  if(.not. allocated(totrefs))  allocate(totrefs(0))
  do ip=1,size(pmh%pc)
    totnver = totnver + size(pmh%pc(ip)%z,2)
    maxdim = max(pmh%pc(ip)%dim, maxdim)
!    do ig = 1, size(pmh%pc(ip)%el,1)
!      maxref = max(maxval(pmh%pc(ip)%el(ig)%ref), maxref)
!    enddo
  enddo

  if(present(z)) then
    if(allocated(z)) deallocate(z)
    allocate(z(maxdim,totnver))
    z = 0
  endif

  maxref = maxref + 1

  print*, 'Writing zone '//trim(string(maxref))//': '//trim(FEDB(check_fe(.true.,1,1,0,0))%desc)

  ini = 1
  write(ini_ch,'(z20)') ini
  write(end_ch,'(z20)') totnver
  write(zone_ch,'(z20)') maxref
  write(dim_ch,'(z20)') maxdim
  write(aux,*) '(10 (',trim(adjustl(zone_ch)), ' ', trim(adjustl(ini_ch)),' ', &
                     & trim(adjustl(end_ch)), ' 1 ',trim(adjustl(dim_ch)),')('
  write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(aux))
  if (ios /= 0) call error('module_write_msh/write_nodes # '//trim(string(ios)))

  cnt = 1
  do ip=1,size(pmh%pc)
    do i=1,size(pmh%pc(ip)%z,2)
      if(present(z)) z(:,cnt) = pmh%pc(ip)%z(:,i)
      write(unit=iu, fmt='(a)', iostat = ios) string(pmh%pc(ip)%z(:,i))
      if (ios /= 0) call error('module_write_msh/write_nodes # '//trim(string(ios)))
      cnt = cnt + 1
    enddo
  enddo

  write(unit=iu, fmt='(a)', iostat = ios) '))'
  if (ios /= 0) call error('module_write_msh/write_nodes # '//trim(string(ios)))


end subroutine


subroutine write_msh_cells(iu, pmh, maxtopdim, maxref)
  integer,           intent(in) :: iu ! file descriptor
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer,           intent(in) :: maxtopdim
  integer,        intent(inout) :: maxref
  integer                       :: ig, ip, i, ios, numcells
  character(len=MAXPATH)        :: ini_ch,end_ch,zone_ch, eltp_ch, aux
  integer, allocatable          :: refs(:)


    numcells = 1
    do ip=1,size(pmh%pc)
      do ig = 1, size(pmh%pc(ip)%el,1)
        if(FEDB(pmh%pc(ip)%el(ig)%type)%tdim == maxtopdim) then
          call sunique(pmh%pc(ip)%el(ig)%ref, refs)
          write(eltp_ch,'(z20)') get_cell_msh_type(pmh%pc(ip)%el(ig)%type)
          do i=1, size(refs,1)
            write(ini_ch,'(z20)') numcells
            numcells = numcells + count(pmh%pc(ip)%el(ig)%ref==refs(i))
            maxref = maxref + 1
            print*, 'Writing zone '//trim(string(maxref))//': '//trim(FEDB(pmh%pc(ip)%el(ig)%type)%desc)
            write(end_ch,'(z20)') numcells-1
            write(zone_ch,'(z20)') maxref
            write(aux,*) '(12 (',trim(adjustl(zone_ch)), ' ', trim(adjustl(ini_ch)),' ', &
                       & trim(adjustl(end_ch)), ' 1 ',trim(adjustl(eltp_ch)),'))'
            write(unit=iu, fmt='(a)', iostat = ios) trim(adjustl(aux))
            if (ios /= 0) call error('module_write_msh/write_nodes # '//trim(string(ios)))
          enddo
  !        call add(alltypes, pmh%pc(ip)%el(ig)%type,size(alltypes), fit=.true.)
        endif
      enddo
    enddo




end subroutine


subroutine build_msh_faces(pmh, maxtopdim)
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer,        intent(in)    :: maxtopdim
  type(msh_faces)               :: faces ! pmh_mesh
  integer                       :: ip, ig, i, j, counter
  integer                       :: tp, lne, lnf, sf, et, ft, fnv, ftd, etd, totfaces, findx
  integer                       :: edge(2,12)
  integer                       :: face(4,6)

  allocate(faces%type(0))
  allocate(faces%zone(0))
  allocate(faces%mm(0))
!  allocate(faces%adcells(0))
  totfaces = 1

  do ip=1, size(pmh%pc,1)
    do ig=1, size(pmh%pc(ip)%el,1)
      if(FEDB(pmh%pc(ip)%el(ig)%type)%tdim == maxtopdim) then
        tp = pmh%pc(ip)%el(ig)%type
        sf = size(faces%type,1)
          counter = 0
        if(maxtopdim == 2) then
          lne = FEDB(tp)%lne
          et = FEDB(pmh%pc(ip)%el(ig)%type)%e_type
          etd = FEDB(FEDB(pmh%pc(ip)%el(ig)%type)%e_type)%tdim
          edge = FEDB(pmh%pc(ip)%el(ig)%type)%edge
          do i=1,pmh%pc(ip)%el(ig)%nel
            do j=1,lne
             findx = face_indx(faces%mm, totfaces, sf+counter, pmh%pc(ip)%el(ig)%mm(edge(:,j),i))
              if(findx < 0) then
                counter = counter + 1
                call add(faces%type, ft, sf+counter, .true.)
                call extend_varlen(faces%mm,sf+counter,.false.)
                faces%mm(sf+counter)%data = pmh%pc(ip)%el(ig)%mm(edge(:,j),i)
                call add(faces%zone, search_face_zone(pmh%pc(ip), pmh%pc(ip)%el(ig)%mm(edge(:,j),i), ftd), sf+counter, .true.)
                call extend_varlen(faces%adcells,sf+counter,.false.)
                call add(faces%adcells(sf+counter)%data, i, 1, .true.)
              else
                call add(faces%adcells(sf+counter)%data, findx, size(faces%adcells(sf+counter)%data,1), .true.)
              endif
            enddo
          enddo
        elseif(maxtopdim == 3) then
          lnf = FEDB(tp)%lnf
          ft = FEDB(pmh%pc(ip)%el(ig)%type)%f_type
          fnv = FEDB(FEDB(pmh%pc(ip)%el(ig)%type)%f_type)%lnv
          ftd = FEDB(FEDB(pmh%pc(ip)%el(ig)%type)%f_type)%tdim
          face = FEDB(pmh%pc(ip)%el(ig)%type)%face
          do i=1,pmh%pc(ip)%el(ig)%nel
            do j=1,lnf
             findx = face_indx(faces%mm, totfaces, sf+counter, pmh%pc(ip)%el(ig)%mm(face(1:fnv,j),i))
              if(findx < 0) then
                counter = counter + 1
                call add(faces%type, ft, sf+counter, .true.)
                call extend_varlen(faces%mm,sf+counter,.false.)
                faces%mm(sf+counter)%data = pmh%pc(ip)%el(ig)%mm(face(1:fnv,j),i)
                call add(faces%zone, search_face_zone(pmh%pc(ip), pmh%pc(ip)%el(ig)%mm(face(1:fnv,j),i), ftd), sf+counter, .true.)
                call extend_varlen(faces%adcells,sf+counter,.false.)
                allocate(faces%adcells(sf+counter)%data(1))
                call add(faces%adcells(sf+counter)%data, i, 1, .true.)
              else
                call add(faces%adcells(sf+counter)%data, findx, size(faces%adcells(sf+counter)%data,1), .true.)
              endif
            enddo
          enddo
        else
          call error('Max topological dimension must be higher than 1')
        endif
        deallocate(pmh%pc(ip)%el(ig)%mm)
        totfaces = sf + counter
      endif
    enddo
  enddo


  call reduce(faces%type, totfaces)
  call reduce(faces%zone, totfaces)
  call reduce_varlen(faces%mm, totfaces)
  call reduce_varlen(faces%adcells, totfaces)
!type msh_faces
!  integer, allocatable :: type(:)
!  integer, allocatable :: zone(:)
!  type(varlen), allocatable :: mm(:)
!  type(varlen), allocatable :: adcells(:)
!end type

end subroutine


!-----------------------------------------------------------------------
! pmh2msh: write a PMH structure into a MSH file
!-----------------------------------------------------------------------
subroutine pmh2msh(iu, pmh, maxref,z)
integer, intent(in)        :: iu ! file descriptor
type(pmh_mesh), intent(in) :: pmh
integer,     intent(inout) :: maxref
real(real64), allocatable, optional, intent(in) :: z(:,:)
integer :: i, ipp, ip, ig, k, j, l, res, max_tdim, nv, nsv, pos, zone, nedges
integer :: numfaces, ini, ios, tempref
integer, allocatable :: piece2save(:), nver_piece(:), nel_piece(:)
integer, allocatable ::  refs(:), new_refs(:), aux_refs(:), ie_zones(:)
character(maxpath) :: str,zone_ch,ini_ch,end_ch,eltype_ch, aux,aux2
logical :: check, ft(2) = [.false.,.true.], tt(2) = [.true.,.true.]
type submmr
  integer              :: n        !number of files of mmr
  integer, allocatable :: tmp(:)   !aux. vector to save an single element (vertices+1)
  integer, allocatable :: mmr(:,:) !matrix mm+ref+2*(adcells) for elements with a dimension smaller than max_tdim (n x (vertices+1))
end type
type(submmr) ::  auxsver(2:4),sver(2:4)         !structure to save elements depending on the number of vertices

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
  call info('(module_pmh/pmh2msh) option -glue not implemented yet')
end if

!calculation of max_tdim: maximal topological dimension
max_tdim = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_pmh/pmh2msh) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  do ig = 1, size(pmh%pc(ip)%el, 1)
    max_tdim = max(max_tdim, FEDB(pmh%pc(ip)%el(ig)%type)%tdim)
  end do
end do

!calculation of sver: mesh faces and edges
if (allocated(nver_piece)) deallocate(nver_piece); !nver_piece(ipp): global numbering for the last vertex of piece #ipp
allocate(nver_piece(0:size(piece2save,1)), stat = res, errmsg = str)
if (res /= 0) call error('(module_pmh/pmh2msh) Unable to allocate variable nver_piece: '//trim(str))
nver_piece(0) = 0
!calculation of sver: mesh faces and edges
if (allocated(nel_piece)) deallocate(nel_piece); !nver_piece(ipp): global numbering for the last vertex of piece #ipp
allocate(nel_piece(0:size(piece2save,1)), stat = res, errmsg = str)
if (res /= 0) call error('(module_pmh/pmh2msh) Unable to allocate variable nel_piece: '//trim(str))
nel_piece(0) = 0
sver(2:4)%n = 0
auxsver(2:4)%n = 0


do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  nver_piece(ipp) = nver_piece(ipp-1) + pmh%pc(ipp)%nver
  nel_piece(ipp) = nel_piece(ipp-1)
  !first round: save in mmr edges and faces (when max_tdim = 3)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if ((FEDB(elg%type)%tdim == 1 .and. max_tdim == 2) .or. (FEDB(elg%type)%tdim == 2 .and. max_tdim == 3)) then
        ! Renumber references
        if(allocated(aux_refs)) deallocate(aux_refs)
        call sunique(pack(elg%ref, elg%ref/=0), aux_refs)
        call ssort(aux_refs)
!        aux_refs = sort(unique(elg%ref))
!        aux_refs = pack(aux_refs,aux_refs/=0) ! Remove reference 0
        if(allocated(aux_refs) .and. size(aux_refs,1)>0) then
          do k = 1, size(aux_refs)
            maxref = maxref + 1
            call set(new_refs, maxref, aux_refs(k))
          enddo
        else
          maxref = maxref+1
        endif
        do k = 1, elg%nel
           nv = FEDB(elg%type)%lnv
          ! If reference is 0
          if(elg%ref(k)==0) then
            tempref = maxref
          else
            tempref = new_refs(elg%ref(k))
          endif
          call set(auxsver(nv)%tmp, [sort(nver_piece(ipp-1) + elg%mm(:,k)), tempref], &
               & [(i, i=1,nv+3)], fit=.true.)
          call insert_sorted(1, auxsver(nv)%mmr, auxsver(nv)%tmp, used=auxsver(nv)%n, fit=ft,pos=pos)
!Nuevo
          sver(nv)%n = auxsver(nv)%n
          call insert(1, sver(nv)%mmr, [nver_piece(ipp-1) + elg%mm(:,k), tempref,0,0],&
            & abs(pos), used=sver(nv)%n, fit=ft)
        enddo
      end if
    end associate
  end do
  !second round: save in mmr edges of faces (when max_tdim = 2) or faces of solids (when max_tdim = 3)
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(elg => pmh%pc(ip)%el(ig)) !elg: current group
      if (max_tdim == 3 .and. FEDB(elg%type)%tdim == 3 ) then !get faces from solids
        print*, 'Building faces from high order elements: '//trim(FEDB(elg%type)%desc)
        ! Renumber references
        if(allocated(aux_refs)) deallocate(aux_refs)
        call sunique(pack(elg%ref, elg%ref/=0), aux_refs)
        call ssort(aux_refs)
!        aux_refs = sort(unique(elg%ref))
!        aux_refs = pack(aux_refs,aux_refs/=0) ! Remove reference 0
        if(allocated(aux_refs) .and. size(aux_refs,1)>0) then
          do k = 1, size(aux_refs)
            maxref = maxref + 1
            call set(new_refs, maxref, aux_refs(k))
          enddo
        else
          maxref = maxref+1
        endif
        if (FEDB(elg%type)%f_type > 0) nsv = FEDB(FEDB(elg%type)%f_type)%lnv !vertices per face (tetr and hexa)
        do k = 1, elg%nel
          do j = 1, FEDB(elg%type)%lnf
            if (FEDB(elg%type)%f_type == 0) nsv = VF_WEDG(j) !vertices per face (wedge)
            ! If reference is 0
            if(elg%ref(k)==0) then
              tempref = maxref
            else
              tempref = new_refs(elg%ref(k))
            endif
            call set(auxsver(nsv)%tmp, [sort(nver_piece(ipp-1) + elg%mm(FEDB(elg%type)%face(1:nsv,j),k)), 0], &
                 & [(i, i=1,nsv+3)], fit=.true.)
            call insert_sorted(1, auxsver(nsv)%mmr, auxsver(nsv)%tmp(1:nsv), used=auxsver(nsv)%n, fit=ft, pos=pos) !ref not saved
            sver(nsv)%n = auxsver(nsv)%n
            if(pos<0) then
!Nuevo
              call insert(1, sver(nsv)%mmr, &
                & [nver_piece(ipp-1) + elg%mm(FEDB(elg%type)%face(nsv:1:-1,j),k),tempref, nel_piece(ipp)+k], &
                & -pos, used=sver(nsv)%n, fit=ft)
!              call set(sver(nsv)%mmr, new_refs(elg%ref(k)), -pos, nsv+1, fit=[.true.,.true.])
!              call set(sver(nsv)%mmr, nel_piece(ipp)+k, -pos, nsv+2, fit=[.true.,.true.])
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checking prints ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print*, '+',nsv,'-', j,k, '-',sver(nsv)%mmr(-pos,nsv+2),'-', trim(string(elg%mm(FEDB(elg%type)%face(nsv:1:-1,j),k)))
            else
              if(sver(nsv)%mmr(pos,nsv+2) == 0) then
                call set(1, sver(nsv)%mmr, [nver_piece(ipp-1) + elg%mm(FEDB(elg%type)%face(nsv:1:-1,j),k)], &
                  & pos, fit=tt)
                call set(sver(nsv)%mmr, nel_piece(ipp)+k, pos, nsv+2, fit=tt)
!print*, '-',nsv,'-', j,k, '-',sver(nsv)%mmr(pos,nsv+2:nsv+3),'-', trim(string(elg%mm(FEDB(elg%type)%face(nsv:1:-1,j),k)))
              else
                call set(sver(nsv)%mmr, nel_piece(ipp)+k, pos, nsv+3, fit=tt)
!print*, '!',nsv,'-', j,k, '-',sver(nsv)%mmr(pos,nsv+2:nsv+3),'-'
              endif
            endif
          end do
        end do
        nel_piece(ipp) = nel_piece(ipp) + elg%nel
      elseif (max_tdim == 2 .and. FEDB(elg%type)%tdim == 2) then !get edges from faces
        print*, 'Building faces from high order elements: '//trim(FEDB(elg%type)%desc)
        ! Renumber references
        if(allocated(aux_refs)) deallocate(aux_refs)
        call sunique(pack(elg%ref, elg%ref/=0), aux_refs)
        call ssort(aux_refs)
!        aux_refs = sort(unique(elg%ref))
!        aux_refs = pack(aux_refs,aux_refs/=0) ! Remove reference 0
        if(allocated(aux_refs) .and. size(aux_refs,1)>0) then
          do k = 1, size(aux_refs)
            maxref = maxref + 1
            call set(new_refs, maxref, aux_refs(k))
          enddo
        else
          maxref = maxref+1
        endif
        if (FEDB(elg%type)%e_type > 0) nsv = FEDB(FEDB(elg%type)%e_type)%lnv !vertices per edge
        do k = 1, elg%nel
          do j = 1, FEDB(elg%type)%lne
            ! If reference is 0
            if(elg%ref(k)==0) then
              tempref = maxref
            else
              tempref = new_refs(elg%ref(k))
            endif
            call set(auxsver(2)%tmp, &
            & [sort(nver_piece(ipp-1) + elg%mm(FEDB(elg%type)%edge(:,j),k)), 0], &
            & [(i, i=1,nsv+3)], .true.)
            call insert_sorted(1, auxsver(nsv)%mmr, auxsver(nsv)%tmp(1:nsv), &
              &used=auxsver(nsv)%n, fit=tt, pos=pos) !ref not saved
            sver(nsv)%n = auxsver(nsv)%n
            if(pos<0) then
!Nuevo
              call insert(1, sver(nsv)%mmr, &
                [nver_piece(ipp-1) + elg%mm(FEDB(elg%type)%edge(:,j),k), tempref, nel_piece(ipp)+k], &
                & -pos, used=sver(nsv)%n, fit=ft) !ref not saved
!              call set(sver(nsv)%mmr, nel_piece(ipp)+k, pos, nsv+2, fit=[.true.,.true.])
            else
              if(sver(nsv)%mmr(pos,nsv+2) == 0) then
                call set(1, sver(nsv)%mmr, [nver_piece(ipp-1) + elg%mm(FEDB(elg%type)%edge(:,j),k)], &
                  & pos, fit=tt)
                call set(sver(nsv)%mmr, nel_piece(ipp)+k, pos, nsv+2, fit=tt)
              else
                call set(sver(nsv)%mmr, nel_piece(ipp)+k, pos, nsv+3, fit=tt)
              endif
            endif
          end do
        end do
        nel_piece(ipp) = nel_piece(ipp) + elg%nel
      end if
    end associate
  end do
end do

numfaces = 0
do i = 2, 4
  sver(i)%n = auxsver(i)%n
  if (allocated(auxsver(i)%mmr)) then
    deallocate(auxsver(i)%mmr)
  endif
  if (allocated(sver(i)%mmr)) then
    call reduce(sver(i)%mmr, sver(i)%n, i+3)
    numfaces = numfaces + sver(i)%n
  endif
enddo


ini = 1
write(ini_ch,'(z20)') ini
write(end_ch,'(z20)') ini+numfaces-1
write(aux,*) '(13 (0 ', trim(adjustl(ini_ch)),' ', &
                      & trim(adjustl(end_ch)), ' 0))'
write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(aux))
if (ios /= 0) call error('module_write_msh/pmh2msh # '//trim(string(ios)))

! Write faces to file
do i = 2, 4
  if (allocated(sver(i)%mmr)) then
    call sunique(unique(sver(i)%mmr(:,i+1)), refs)
    call ssort(refs)
!    refs = sort(unique(sver(i)%mmr(:,i+1)))
    write(eltype_ch,'(z20)') i
    do j=1, size(refs,1)
      check = .false.
      numfaces = count( sver(i)%mmr(:,i+1)==refs(j) )
      if(refs(j) == 0) then
        maxref = maxref + 1
        zone = maxref
      else
        zone = refs(j)
      endif
      write(ini_ch,'(z20)') ini
      write(end_ch,'(z20)') ini+numfaces-1
      write(zone_ch,'(z20)') zone
      if (ios /= 0) call error('module_write_msh/pmh2msh # '//trim(string(ios)))
      if(i==2) then; nedges = 1; else; nedges = i; endif
      print*, 'Writing zone '//trim(string(zone))//': '//trim(FEDB(check_fe(.true.,i,i,nedges,0))%desc)
      do k=1, size(sver(i)%mmr,1)
        if(sver(i)%mmr(k,i+1)==refs(j)) then
          ! interior/exterior groups check
          if(.not. check) then
            if(size(unique(sver(i)%mmr(:,i+3)),1) /= 1 .and. sver(i)%mmr(k,i+3) == 0) then
              ! exterior zone: wall
              call set(ie_zones, 3, j)
            else
              ! interior zone
              call set(ie_zones, 2, j)
            endif
            ! writes face group header
            write(aux,*) '(13 (',trim(adjustl(zone_ch)),' ', trim(adjustl(ini_ch)),' ', &
                         & trim(adjustl(end_ch)), ' ', trim(string(ie_zones(j))),' ', &
                         & trim(adjustl(eltype_ch)),')('
            write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(aux))
            check = .true.
          endif
          aux = ''
!          if(present(z)) call face_positive_jacobian(sver(i)%mmr(k,:i), z)
          do l=1,i
            write(aux2,'(z20)') sver(i)%mmr(k,l)
            write(aux,'(a)') trim(adjustl(aux))//' '//trim(adjustl(aux2))//' '
          enddo
          write(aux2,'(z20)') sver(i)%mmr(k,i+2)
          write(aux,'(a)') trim(adjustl(aux))//' '//trim(adjustl(aux2))//' '
          write(aux2,'(z20)') sver(i)%mmr(k,i+3)
          write(aux,'(a)') trim(adjustl(aux))//' '//trim(adjustl(aux2))//' '
          write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(aux))
          if (ios /= 0) call error('module_write_msh/pmh2msh # '//trim(string(ios)))
        endif
      enddo
      ini = ini+numfaces
      write(unit=iu, fmt='(a)', iostat = ios)  '))'
      if (ios /= 0) call error('module_write_msh/pmh2msh # '//trim(string(ios)))
    enddo
  endif
end do

! Write zones (45
do j=1,size(refs,1)
  if(ie_zones(j)==2) then
    write(aux,'(a)') '(45 ('//trim(string(refs(j)))//' interior interior'//trim(string(j))//')())'
    write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(aux))
    if (ios /= 0) call error('module_write_msh/pmh2msh # '//trim(string(ios)))
  elseif(ie_zones(j)==3) then
    write(aux,'(a)') '(45 ('//trim(string(refs(j)))//' wall wall'//trim(string(j))//')())'
    write(unit=iu, fmt='(a)', iostat = ios)  trim(adjustl(aux))
    if (ios /= 0) call error('module_write_msh/pmh2msh # '//trim(string(ios)))
  else
    call info('Unknown zone type: '//trim(string(string(refs(j)))))
  endif
enddo

end subroutine


function face_indx(mm, ini, s, el_mm) result(res)
  type(varlen),   intent(in) :: mm(:)
  integer,        intent(in) :: ini      ! indes where to start to search
  integer,        intent(in) :: s        ! number of elements in mm
  integer,        intent(in) :: el_mm(:) ! element connectivity
  integer                    :: res
  integer                    :: i

  res = -1

  do i=ini, s
    if(all(bubble(mm(i)%data) .eq. bubble(el_mm))) then
      res = i
      return
    endif
  enddo

end function


function search_face_zone(pc, el_mm, tdim) result(zone)
  type(piece),    intent(in) :: pc      ! pmh piece
  integer,        intent(in) :: el_mm(:) ! element connectivity
  integer,        intent(in) :: tdim
  integer                    :: zone
  integer                    :: ig, i

  zone = 0
  do ig=1, size(pc%el,1)
    if(FEDB(pc%el(ig)%type)%tdim == tdim) then
      do i=1, pc%el(ig)%nel
        if(all(bubble(pc%el(ig)%mm(:,i)) .eq. bubble(el_mm))) then
            zone = pc%el(ig)%ref(i)
            return
        endif
      enddo
    endif
  enddo

end function


end module
