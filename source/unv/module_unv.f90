module module_unv_fcnv
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
use basicmod, only: real64, maxpath, error, info, string, int, word_count, lcase, sfind, sunique, is_arg, get_post_arg
use module_pmh_fcnv
use module_fe_database_pmh_fcnv, only: FEDB, check_fe
use module_manage_unv_fcnv
use module_mesh_unv_fcnv
implicit none

contains

!-----------------------------------------------------------------------
! load_unv: read a UNV file
!-----------------------------------------------------------------------
subroutine load_unv(unvfile, pmh, padval, infield, ca_opt)
character(len=*), intent(in) :: unvfile
type(pmh_mesh),intent(inout) :: pmh
real(real64),     intent(in) :: padval
character(len=*), allocatable, intent(in) :: infield(:)
logical,          intent(in) :: ca_opt
type(unv)                    :: u

!inital settings
!call report_option('level', 'stdout')
!process universal file
call open_unv(u, unvfile)
call read_unv(u, pmh, padval, infield, ca_opt)
call build_vertices(pmh)
end subroutine

!-----------------------------------------------------------------------
! save_unv: save a PMH structure into a UNV file
!-----------------------------------------------------------------------
subroutine save_unv(filename, iu, pmh, infield, outfield, ca_opt)
character(*),   intent(in) :: filename !mesh filename
integer,        intent(in) :: iu       !file unit
type(pmh_mesh), intent(in) :: pmh      !pmh structure
character(len=*), allocatable, intent(in) :: infield(:) !List of input field names
character(len=*), allocatable, intent(in) :: outfield(:) !List of output field names
logical,        intent(in) :: ca_opt

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
call info('Writing coordinates ...')
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
  if (any_P2) then !save node coordinates when there is some Lagrange P2 elements
    do j = 1, pmh%pc(ip)%nnod
      write(iu,'(4I10)') prev_coord + j, 0, 0, 11
      write(iu,'(1P3E25.16)') reshape(znod(:,j), [3], [0._real64,0._real64,0._real64])
    end do
  else !save vertex coordinates when there is not any Lagrange P2 element
    do j = 1, pmh%pc(ip)%nver
      write(iu,'(4I10)') prev_coord + j, 0, 0, 11
      write(iu,'(1P3E25.16)') reshape(pmh%pc(ip)%z(:,j), [3], [0._real64,0._real64,0._real64])
    end do
  end if
  !update the number of previous nodes/vertices
  if (any_P2) then; prev_coord = prev_coord + pmh%pc(ip)%nnod
  else;             prev_coord = prev_coord + pmh%pc(ip)%nver
  end if
end do
write(iu,'(I6)') -1

!save elements
call info('Writing conectivities ...')
call FE_DB_init() !initialize FE database
write(iu,'(I6)') -1
write(iu,'(I6)') 2412
prev_coord = 0; prev_nel = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
      call info('  Element type: '//trim(FEDB(tp)%desc))
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
call info('Writing references ...')
write(iu,'(I6)') -1
write(iu,'(I6)') 2467
prev_coord = 0; prev_nel = 0
do ipp = 1, size(piece2save,1)
  ip = piece2save(ipp)
  do ig = 1, size(pmh%pc(ip)%el,1)
    associate(elg => pmh%pc(ip)%el(ig), tp => pmh%pc(ip)%el(ig)%type)
      call sunique(elg%ref, unique_ref)
      do ir = 1, size(unique_ref,1)
        call info('  piece_'//trim(string(ip))//'_group_'//trim(string(ig))//'_ref_'//trim(string(unique_ref(ir))))
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

!save fields
call write_unv_fields(iu, pmh, piece2save, infield, outfield, ca_opt)
close(iu)

end subroutine

!-----------------------------------------------------------------------
! save_unv: save a PMH field into a UNV file
!-----------------------------------------------------------------------
subroutine write_unv_fields(iu, pmh, piece2save, infield, outfield, ca_opt, dataset, nparam)
integer,         intent(in) :: iu       !file unit
type(pmh_mesh),  intent(in) :: pmh      !pmh mesh structure
integer,         intent(in) :: piece2save(:)
character(len=*),allocatable, intent(in) :: infield(:) ! List of input field names
character(len=*),allocatable, intent(in) :: outfield(:) ! List of output field names
logical,         intent(in) :: ca_opt   !Code Aster option
integer,optional,intent(in) :: dataset
integer,optional,intent(in) :: nparam
integer                     :: prev_coord
integer                     :: prev_nel
integer                     :: i, j, k, l, m
integer                     :: counter, ds, ncomp, dc, np
character(len=maxpath)      :: fieldname

if(present(dataset)) then; ds = dataset
else;                      ds = 2414
endif
if(present(nparam)) then; np = nparam
else;                     np = 1
endif

prev_coord = 0; prev_nel = 0; counter = 0
do i = 1,size(piece2save,1)
  associate(pc => pmh%pc(piece2save(i)))
    ! Check if exits imposed field names
    if(allocated(outfield)) then
      if(allocated(infield)) then
        if(size(infield,1) /= size(outfield,1)) &
          & call error('Number of input and output field names must agree.')
      else
        if(size(outfield,1) /= get_piece_num_fields(pc)) &
          & call error('Number of output field names must agree with number of fields.')
      endif
    endif
    if(ca_opt) ds = 55  ! Code Aster option forces to write dataset 55 for node fields
    if(allocated(pc%fi) .and. (ds==2414 .or. ds==55)) then
      ! Node fields
      do j = 1,size(pc%fi,1)
        counter = counter + 1
        fieldname =  trim(adjustl(pc%fi(j)%name))
        ! Check if exits imposed field names
        if(allocated(outfield)) then
          if(allocated(infield)) then
            do k=1,size(infield,1)
              if(trim(adjustl(pc%fi(j)%name))==trim(adjustl(infield(k)))) &
                & fieldname = trim(adjustl(outfield(k)))
            enddo
          else
            fieldname = trim(adjustl(outfield(k)))
          endif
        endif
        call info('Writing node field: '//trim(adjustl(fieldname)))
        ncomp = size(pc%fi(j)%val,1)
        if(ncomp == 1) then;     dc = 1 ! Scalar
        elseif(ncomp == 3) then; dc = 2 ! 3 Comp. vector
        elseif(ncomp == 6) then; dc = 4 ! 6 Comp. simmetric global tensor
        elseif(ncomp == 9) then; dc = 5 ! 9 Comp. general tensor
        else;                    dc = 0
        endif
        do np=1,size(pc%fi(j)%param,1)
          write(iu,'(1I6)') -1
          write(iu,'(1I6)') ds                          ! Dataset
          if(ds == 2414) then
            write(iu,'(1I10)') counter                  ! Record1: Analysis dataset label (2414)
            write(iu,*) trim(adjustl(fieldname))        ! Record2: Analysis dataset name (2414)
            write(iu,'(1I10)') 1                        ! Record3: Dataset location. 1:Data at nodes, 2:Data on elements (2414)
          endif
          write(iu,*) trim(adjustl(fieldname))          ! Record4: ID line 1 (2414)
          write(iu,*) 'Double precision floating point' ! Record5: ID line 2 (2414)
          write(iu,*) 'NONE'                            ! Record6: ID line 3 (2414)
          write(iu,*) 'NONE'                            ! Record7: ID line 4 (2414)
          write(iu,*) 'NONE'                            ! Record8: ID line 5 (2414)
          ! Model type, Analysis type, Data characteristic, Result type, Data type, Number of data values for data component
          ! Record9: Unknown, Unknown, dc, User defined, Double precision, ncomp (2414)
          if(ca_opt .and. is_ca_field_type(pc%fi(j)%name)) then
            write(iu,'(6I10)') get_ca_field_record9(pc%fi(j)%name, (/0,0,dc,1000+counter,4,ncomp/))
          else
            write(iu,'(6I10)') 0,0,dc,1000+counter,4,ncomp
          endif
          if(ds == 2414) then
            !Design set ID, Iteration number, Solution set ID, Boundary condition, Load set, Mode number, Time step number,
            !Frequency number, ...
            write(iu,'(8I10)') (/j,np-1,j,0,0,np-1,np-1,np-1/) ! Record10: Integer analysis type speciic data (1-8) (2414)
            write(iu,'(8I10)') (/(0,m=1,2)/)                   ! Record11: Integer analysis type speciic data (9-10) (2414)
            !EigenValue, Modal Mass, ...
            write(iu,'(6E13.5)') &                             ! Record12: Real analysis type specific data (1-6) (2414)
              & (/pc%fi(j)%param(np),pc%fi(j)%param(np),(0._real64,m=1,4)/)
            write(iu,'(6E13.5)') (/(0._real64,m=1,6)/)         ! Record13: Real analysis type specific data (7-12) (2414)
          else
            write(iu,'(8I10)') ncomp,ncomp,j,np-1,0,0,0,0      ! Record7: Integer analysis type speciic data (1-8) (55)
            write(iu,'(6E13.5)') &                             ! Record8: Real analysis type specific data (9-10) (55)
              & (/pc%fi(j)%param(np),(0._real64,m=1,5)/)
          endif
          do k = 1, pc%nnod
            write(iu,'(1I10)') prev_coord+k                    ! Record14: Node Number (2414)
            write(iu,'(6E13.5)') &                             ! Record15: Data at this node (2414)
              & reshape(pc%fi(j)%val(:,k,np), [ncomp], [(0._real64,m=1,ncomp)])
          enddo
          write(iu,'(I6)') -1
        enddo
        prev_coord = prev_coord+pc%nnod
      enddo
    endif
    do j=1, size(pc%el,1)
      if(FEDB(pc%el(j)%type)%tdim == 0) cycle
      if(ca_opt) ds = 57  ! Code Aster option forces to write dataset 57 for node fields
      if(allocated(pc%el(j)%fi) .and. (ds==2414 .or. ds==57)) then
        do k=1, size(pc%el(j)%fi,1)
          counter = counter + 1
          fieldname =  trim(adjustl(pc%el(j)%fi(k)%name))
          ! Check if exits imposed field names
          if(allocated(outfield)) then
            if(allocated(infield)) then
              do l=1,size(infield,1)
                if(trim(adjustl(pc%el(j)%fi(k)%name))==trim(adjustl(infield(l)))) &
                  & fieldname = trim(adjustl(outfield(l)))
              enddo
            else
              fieldname = trim(adjustl(outfield(counter)))
            endif
          endif
          call info('Writing cell field: '//trim(adjustl(fieldname)))
          ncomp = size(pc%el(j)%fi(k)%val,1)

          if(ncomp == 1) then;                                dc = 1 ! Scalar
          elseif(ncomp == 2 .or. ncomp == 3) then; ncomp = 3; dc = 2 ! 3 Comp. vector
          elseif(ncomp == 6) then;                            dc = 4 ! 6 Comp. simmetric global tensor
          elseif(ncomp == 9) then;                            dc = 5 ! 9 Comp. general tensor
          else;                                               dc = 0
          endif
          do np=1,size(pc%el(j)%fi(k)%param,1)
            write(iu,'(1I6)') -1
            write(iu,'(1I6)') ds                          ! Dataset
            if(ds == 2414) then
              write(iu,'(1I10)') counter                  ! Record1: Analysis dataset label (2414)
              write(iu,*) trim(adjustl(fieldname))        ! Record2: Analysis dataset name (2414)
              write(iu,'(1I10)') 2                        ! Record3: Dataset location. 1:Data at nodes, 2:Data on elements (2414)
            endif
            write(iu,*) trim(adjustl(fieldname))          ! Record4: ID line 1 (2414)
            write(iu,*) 'Double precision floating point' ! Record5: ID line 2 (2414)
            write(iu,*) 'NONE'                            ! Record6: ID line 3 (2414)
            write(iu,*) 'NONE'                            ! Record7: ID line 4 (2414)
            write(iu,*) 'NONE'                            ! Record8: ID line 5 (2414)
            ! Model type, Analysis type, Data characteristic, Result type, Data type, Number of data values for data component
            ! Record9: Unknown, Unknown, dc, User defined, Double precision, ncomp ( Default, 2414)
            if((ca_opt) .and. &
              & is_ca_field_type(fieldname)) then
              write(iu,'(6I10)') &
                & get_ca_field_record9(fieldname, (/0,0,dc,1000+counter,4,ncomp/))
            else
              write(iu,'(6I10)') 0,0,dc,1000+counter,4,ncomp
            endif
            if(ds == 2414) then
              write(iu,'(8I10)') (/(0,m=1,8)/)                 ! Record10: Integer analysis type speciic data (1-8) (2414)
              write(iu,'(8I10)') (/(0,m=1,2)/)                 ! Record11: Integer analysis type speciic data (9-10) (2414)
              write(iu,'(6E13.5)') (/(0._real64,m=1,6)/)       ! Record12: Real analysis type speciic data (1-6) (2414)
              write(iu,'(6E13.5)') (/(0._real64,m=1,6)/)       ! Record13: Real analysis type speciic data (7-12) (2414)
            else
              write(iu,'(8I10)') ncomp,ncomp,k,np-1,0,0,0,0    ! Record7: Integer analysis type speciic data (9-10) (57)
              write(iu,'(6E13.5)') &                           ! Record8: Real analysis type specific data (9-10) (57)
              & (/pc%el(j)%fi(k)%param(np),(0._real64,m=1,5)/)
            endif
            do l = 1, pc%el(j)%nel
              if(ds == 2414) then                              ! Record14: Element Number, Number of data values (2414)
                write(iu,'(2I10)') prev_nel+l,ncomp
              else
                write(iu,'(2I10)') prev_nel+l
              endif
              write(iu,'(6E13.5)') &                           ! Record15: Data at this element
                & reshape(pc%el(j)%fi(k)%val(:,l,np), [ncomp], [(0._real64,m=1,ncomp)])
            enddo
            write(iu,'(I6)') -1
          enddo
        enddo
      endif
      prev_nel = prev_nel+pc%el(j)%nel
    enddo
  end associate
enddo
end subroutine

end module
