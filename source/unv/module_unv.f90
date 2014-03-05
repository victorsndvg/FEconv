module module_unv
!-----------------------------------------------------------------------
! Module to manage UNV (I-Deas universal) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 27/05/2013
!
! PUBLIC PROCEDURES:
!   load_unv: loads a mesh from a UNV format file
!   save_unv: save a PMH structure into a UNV file
!
!   NOTES: 1) mesh must be composed of a single FE type
!          2) planar meshes must lay in the XY plane
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_convers, only: string, int
use module_alloc, only: sfind
use module_set, only: sunique
use module_args, only: is_arg, get_post_arg
use module_pmh, only: pmh_mesh, build_node_coordinates
use module_fe_database_pmh, only: FEDB, check_fe
use module_manage_unv
use module_mesh
implicit none

contains

!-----------------------------------------------------------------------
! load_unv: read a UNV file
!-----------------------------------------------------------------------
subroutine load_unv(unvfile, nel, nnod, nver, DI_, LNN, LNV, LNE, LNF, nn, mm, nrc, nra, nrv, z, nsd, is_opt)
character(len=*),         intent(in)  :: unvfile
integer,                  intent(out) :: nel, nnod, nver, DI_, LNN, LNV, LNE, LNF
integer, allocatable,     intent(out) :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64),allocatable, intent(out) :: z(:,:)
logical, intent(in) :: is_opt
type(unv):: u
type(mfm_mesh), allocatable :: m(:)
integer :: d, DIM

DIM = 3 !it can change after reading dataset 2412
!inital settings
call report_option('level', 'stdout')
!allocate mesh(es)
allocate(m(DIM))
forall (d = 1:DIM) m(d)%DIM = d
!process universal file
call open_unv(u, unvfile)
call read_unv(u, m, DIM, is_opt)

!return msh variables
nel  = m(DIM)%nl
nnod = m(DIM)%nd
nver = m(DIM)%nv
DI_  = m(DIM)%DIM
LNN  = m(DIM)%LNN
LNV  = m(DIM)%LNV
LNE  = m(DIM)%LNE
LNF  = m(DIM)%LNF
if (allocated(m(DIM)%id)) then; allocate( nn(LNN,nel));  nn(1:LNN,1:nel)  = m(DIM)%id; end if
if (allocated(m(DIM)%iv)) then; allocate( mm(LNV,nel));  mm(1:LNV,1:nel)  = m(DIM)%iv; end if
if (allocated(m(DIM)%rf)) then; allocate(nrc(LNF,nel)); nrc(1:LNF,1:nel)  = m(DIM)%rf; end if
if (allocated(m(DIM)%re)) then; allocate(nra(LNE,nel)); nra(1:LNE,1:nel)  = m(DIM)%re; end if
if (allocated(m(DIM)%rv)) then; allocate(nrv(LNV,nel)); nrv(1:LNV,1:nel)  = m(DIM)%rv; end if
if (allocated(m(DIM)%xv)) then; allocate(  z(3,nver));  z=0.0_real64; z(1:DI_,1:nver) = m(DIM)%xv; DI_=3; end if
if (allocated(m(DIM)%rl)) then; allocate(    nsd(nel));       nsd(1:nel)  = m(DIM)%rl; end if

!deallocate mesh
do d = 1, DIM
  call dealloc_mesh(m(d))
enddo
deallocate(m)
end subroutine

!-----------------------------------------------------------------------
! save_unv: save a PMH structure into a UNV file
!-----------------------------------------------------------------------
subroutine save_unv(filename, iu, pmh)
character(*),   intent(in) :: filename !mesh filename
integer,        intent(in) :: iu       !file unit
type(pmh_mesh), intent(in) :: pmh      !pmh structure

integer :: i, ipp, ip, ig, ios, prev_nel, prev_coord, j, k, idesc, l, ir
integer, allocatable :: piece2save(:), unique_ref(:), k_ref(:)
real(real64), allocatable :: znod(:,:)
character(maxpath) :: str
logical :: any_P2, any_RT_ND, all_are_P1

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
  call info('(module_pmh/pmh2mfm) option -glue not implemented yet')
end if

!open file
open (unit = iu, file = filename, form = 'formatted', position = 'rewind', iostat = ios)
if (ios /= 0) call error('(module_unv/save_unv) unable to open file; error number '//trim(string(ios)))

!save coordinates
write(iu,'(I6)') -1
write(iu,'(I6)') 2411
prev_coord = 0; prev_nel = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  if (1 > ip .or. ip > size(pmh%pc, 1)) call error('(module_unv/save_unv) requested piece '//trim(string(ip))//&
  &' does not exist in the mesh')
  !check whether any element group is Lagrange P2, Raviart-Thomas (edge) or Nedelec (face)
  any_P2    = .false.
  any_RT_ND = .false.
  do ig = 1, size(pmh%pc(ip)%el, 1)
    associate(tp => pmh%pc(ip)%el(ig)%type)
      if (.not. FEDB(tp)%nver_eq_nnod .and. (FEDB(tp)%lnn == FEDB(tp)%lnv + FEDB(tp)%lne)) any_P2 = .true.
      if (tp == check_fe(.false.,  3, 3,  3, 0) .or. &
          tp == check_fe(.false.,  6, 4,  6, 4) .or. &
          tp == check_fe(.false.,  4, 4,  6, 4)) any_RT_ND = .true.
    end associate
  end do
  if (any_P2 .and. any_RT_ND) call error('(module_unv/save_unv) piece '//trim(string(ip))//' has groups with '//&
  &'Lagrange P2 and edge or face elements mixed; unable to save mesh in UNV format')
  call build_node_coordinates(pmh%pc(ip), ip, all_are_P1, znod)
  if (any_P2) then !save node coordinates (if there is any Lagrange P2 element)
    do j = 1, pmh%pc(ip)%nnod
      write(iu,'(4I10)') prev_coord + j, 0, 0, 11
      write(iu,'(1P3D25.16)') reshape(znod(:,j), [3], [0._real64,0._real64,0._real64])
    end do
  else !save vertex coordinates (if there is not any Lagrange P2 element)
    do j = 1, pmh%pc(ip)%nver
      write(iu,'(4I10)') prev_coord + j, 0, 0, 11
      write(iu,'(1P3D25.16)') reshape(pmh%pc(ip)%z(:,j), [3], [0._real64,0._real64,0._real64])
    end do
  end if
  !update the number of previous nodes/vertices
  if (any_P2) then; prev_coord = prev_coord + pmh%pc(ip)%nnod
  else;             prev_coord = prev_coord + pmh%pc(ip)%nver
  end if
end do
write(iu,'(I6)') -1

!save elements
call FE_DB_init() !initialize FE database
write(iu,'(I6)') -1
write(iu,'(I6)') 2412
prev_coord = 0; prev_nel = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
      if (FEDB(tp)%tdim > 1) then !non-beam elements
        idesc = check_unv_fe(FEDB(tp)%tdim, FEDB(tp)%lnn, FEDB(tp)%lnv, FEDB(tp)%lne, FEDB(tp)%lnf)
        if (idesc == 0) call error('(module_unv/save_unv) unable to find a UNV equivalente to element type '//&
        &trim(FEDB(tp)%desc)//'; piece '//trim(string(ip))//', group '//trim(string(ig)))
        if (FEDB(tp)%lnn == FEDB(tp)%lnv + FEDB(tp)%lne) then !Lagrange P2 elements      
          do k = 1, elg%nel
            write(iu,'(6I10)') prev_nel + k, idesc, 2, 1, 7, FEDB(tp)%lnn
            do l = 0, (FEDB(tp)%lnn-1)/8
              write(iu,'(8I10)') (prev_coord + elg%nn(FE_DB(idesc)%nn_order(i),k), i = 1+8*l, min(8*(l+1), FEDB(tp)%lnn))
            end do 
          end do
        else !non-Lagrange P2 elements      
          do k = 1, elg%nel
            write(iu,'(6I10)') prev_nel + k, idesc, 2, 1, 7, FEDB(tp)%lnv
            do l = 0, (FEDB(tp)%lnv-1)/8
              !ATTENTION: if we have P1 and P2 mixed, znod was written above; we assume here that mm has the same global numbering 
              !than nn in vertices
              write(iu,'(8I10)') (prev_coord + elg%mm(i,k), i = 1+8*l, min(8*(l+1), FEDB(tp)%lnv)) 
            end do 
          end do
        end if
        prev_nel = prev_nel + elg%nel
      elseif (FEDB(tp)%tdim == 1) then !beam elements
        idesc = check_unv_fe(FEDB(tp)%tdim, FEDB(tp)%lnn, FEDB(tp)%lnv, FEDB(tp)%lne, FEDB(tp)%lnf)
        if (idesc == 0) call error('(module_unv/save_unv) find a UNV correspondence to element type '//&
        &trim(FEDB(tp)%desc)//'; piece '//trim(string(ip))//', group '//trim(string(ig)))
        if (FEDB(tp)%lnn == FEDB(tp)%lnv + FEDB(tp)%lne) then !Lagrange P2 elements      
          do k = 1, elg%nel
            write(iu,'(6I10)') prev_nel + k, idesc, 2, 1, 7, FEDB(tp)%lnn
            write(iu,'(3I10)') 0, 0, 0 
            write(iu,'(8I10)') (prev_coord + elg%nn(FE_DB(idesc)%nn_order(i),k), i = 1, FEDB(tp)%lnn)
          end do
        else !non-Lagrange P2 elements      
          do k = 1, elg%nel
            write(iu,'(6I10)') prev_nel + k, idesc, 2, 1, 7, FEDB(tp)%lnv
            write(iu,'(3I10)') 0, 0, 0 
            !ATTENTION: if we have P1 and P2 mixed, znod was written above; we assume here that mm has the same global numbering 
            !than nn in vertices
            write(iu,'(8I10)') prev_coord + elg%mm(:,k)
          end do
        end if
        prev_nel = prev_nel + elg%nel
      end if !for tdim = 0 (vertices), do not update prev_nel
    end associate
  end do
  !update the number of previous nodes/vertices
  if (any_P2) then; prev_coord = prev_coord + pmh%pc(ip)%nnod
  else;             prev_coord = prev_coord + pmh%pc(ip)%nver
  end if
end do
write(iu,'(I6)') -1

!save groups
write(iu,'(I6)') -1
write(iu,'(I6)') 2467
prev_coord = 0; prev_nel = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
      call sunique(elg%ref, unique_ref)
      do ir = 1, size(unique_ref,1)
        call sfind(elg%ref, unique_ref(ir), k_ref)
        write(iu,'(8I10)') unique_ref(ir), 0, 0, 0, 0, 0, 0, size(k_ref, 1) 
        write(iu,'(A)') 'piece_'//trim(string(ip))//'_group_'//trim(string(ig))//'_ref_'//trim(string(unique_ref(ir)))
        if (FEDB(tp)%tdim > 0) then !non-vertex elements
          write(iu,'(8i10)') (8, prev_nel + k_ref(i), 0, 0, i = 1, size(k_ref,1))
        else !vertex elements
          write(iu,'(8i10)') (7, prev_coord + elg%mm(1,k_ref(i)), 0, 0, i = 1, size(k_ref,1))
        end if
      end do
      if (FEDB(tp)%tdim > 0) prev_nel = prev_nel + elg%nel
    end associate
  end do
  !update the number of previous nodes/vertices
  if (any_P2) then; prev_coord = prev_coord + pmh%pc(ip)%nnod
  else;             prev_coord = prev_coord + pmh%pc(ip)%nver
  end if
end do
write(iu,'(I6)') -1
close(iu)
end subroutine

end module
