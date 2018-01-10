module module_mff_fcnv
!-----------------------------------------------------------------------
! Module to manage MFF (Modulef Formatted Field) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Víctor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/05/2014
!
! PUBLIC PROCEDURES:
!   load_mff: loads a mesh from a MFM format file
!   save_mff: saves a mesh in a MFM format file
!-----------------------------------------------------------------------
use basicmod, only: real64, error, string, replace, get_unit, feed, empty
use module_pmh_fcnv

implicit none

!type field
!  character(maxpath)        :: name
!  real(real64), allocatable :: param(:)   !nshot
!  real(real64), allocatable :: val(:,:,:) !ncomp x nnod x nshot
!end type
!private :: field

contains

subroutine load_mff(pmh, filenames, fieldnames, param)
  character(len=*), allocatable, intent(in) :: filenames(:) !fields file names
  type(pmh_mesh),             intent(inout) :: pmh
  character(*), allocatable                 :: fieldnames(:) !field names
  real(real64), optional,        intent(in) :: param
  character(len=maxpath)                    :: filename, fieldname
  integer                                   :: iu, ios, i, j, k, idx
  integer                                   :: ncomp, totcomp, maxtdim, aux
  real(real64), allocatable                 :: fielddata(:)
  type(field), allocatable                  :: tempfields(:)


  if(.not. allocated(filenames)) call error('load_mff, filenames is not allocated.')
  if(.not. allocated(fieldnames)) then
    call info('load_mff, fieldnames is not allocated; use filenames instead.')
    call alloc(fieldnames, size(filenames))
    fieldnames = filenames
  end if
  if(size(filenames,1) /= size(fieldnames,1)) call error('load_mff, dimension of filenames and fieldnames must agree.')
  do j = 1, size(filenames,1)
    filename = trim(filenames(j))
    fieldname = trim(fieldnames(j))
    !open file
    iu = get_unit()
    open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
    if (ios /= 0) call error('load_mff/open, #'//trim(string(ios)))
    !read number of components, totcomp
    read(unit=iu, fmt=*, iostat=ios) totcomp
    backspace(unit=iu, iostat=ios)
    if (ios /= 0) call error('load_mff/backspace, #'//trim(string(ios)))
    if(allocated(fielddata)) deallocate(fielddata)
    allocate(fielddata(totcomp))
    !read components, fielddata
    read(unit=iu, fmt=*, iostat=ios) aux, (fielddata(i),  i=1,totcomp)
    close(iu)

    maxtdim = 0
    if(size(pmh%pc,1) == 1) then
      !field over nodes. 1,2 or 3 components allowed
      if (mod(totcomp, size(pmh%pc(1)%z,2)) == 0 .and. (totcomp/size(pmh%pc(1)%z,2)<= 3)) then
        ncomp = totcomp/size(pmh%pc(1)%z,2)
        if(.not. allocated(pmh%pc(1)%fi)) then
          !allocate fi for the first time
          allocate(pmh%pc(1)%fi(1))
        else
          !increase the size of fi in one
          if(allocated(tempfields)) deallocate(tempfields)
          allocate(tempfields(size(pmh%pc(1)%fi,1)+1))
          tempfields(1:size(pmh%pc(1)%fi,1)) = pmh%pc(1)%fi(:)
          call move_alloc(from=tempfields, to=pmh%pc(1)%fi)
        endif
        !current field, idx
        idx = size(pmh%pc(1)%fi,1)
        call info('Reading node field "'//trim(adjustl(fieldname))//'" from: '//trim(adjustl(filename)))
        associate(fix => pmh%pc(1)%fi(idx))
          fix%name = trim(fieldname)
          if(allocated(fix%param)) deallocate(fix%param)
          allocate(fix%param(1))
          if(present(param)) then
            fix%param(1) = param
          else
            fix%param(1) = 0._real64
          endif
          if(allocated(fix%val)) deallocate(fix%val)
          allocate(fix%val(ncomp, totcomp/ncomp,1))
          do i=1, totcomp/ncomp
            fix%val(:,i,1) = fielddata((i-1)*ncomp+1:i*ncomp)
          enddo
        end associate
      else
        !possible celldata, first detect in which elgroup is
        associate(pc1 => pmh%pc(1))
          do i = 1, size(pc1%el,1)
            maxtdim = max(FEDB(pc1%el(i)%type)%tdim, maxtdim)
          enddo
          do i = 1, size(pc1%el,1)
            if(maxtdim == FEDB(pc1%el(i)%type)%tdim) then
              if(mod(totcomp, pc1%el(i)%nel) == 0 .and. totcomp/pc1%el(i)%nel<= 3) then
                ncomp = totcomp/pc1%el(i)%nel
                if(.not. allocated(pc1%el(i)%fi)) then
                  allocate(pc1%el(i)%fi(1))
                else
                  if(allocated(tempfields)) deallocate(tempfields)
                  allocate(tempfields(size(pc1%el(i)%fi,1)+1))
                  tempfields(1:size(pc1%fi,1)) = pc1%el(i)%fi(:)
                  call move_alloc(from=tempfields, to=pc1%el(i)%fi)
                  deallocate(tempfields)
                endif
                idx = size(pc1%el(i)%fi,1)
                call info('Reading cell field "'//trim(adjustl(fieldname))//'" from: '//trim(adjustl(filename)))
                pc1%el(i)%fi(idx)%name = trim(fieldname)
                if(allocated(pc1%el(i)%fi(idx)%param)) deallocate(pc1%el(i)%fi(idx)%param)
                allocate(pc1%el(i)%fi(idx)%param(1))
                if(present(param)) then
                  pc1%el(i)%fi(idx)%param(1) = param
                else
                  pc1%el(i)%fi(idx)%param(1) = 0._real64
                endif
                if(allocated(pc1%el(i)%fi(idx)%val)) deallocate(pc1%el(i)%fi(idx)%val)
                allocate(pc1%el(i)%fi(idx)%val(ncomp, totcomp/ncomp,1))
                do k=1, pc1%el(i)%nel
                  pc1%el(i)%fi(idx)%val(:,k,1) = fielddata((k-1)*ncomp+1:k*ncomp)
                enddo
                exit
              endif
            endif
          enddo
        end associate
      endif
    endif
  enddo
  print '(a)', 'Done!'
end subroutine

subroutine save_mff(pmh, outfield, path, param)
  type(pmh_mesh),            intent(inout) :: pmh      !PMH mesh
  character(*), allocatable, intent(in) :: outfield(:) !Out field file names
  character(*),              intent(in) :: path        !Path
  real(real64), optional,    intent(in) :: param
  character(len=maxpath)                :: filename    !File names
  integer                               :: i,j,k,l,m,pi,mtdim
  integer                               :: iu, ios, fidx

  if(.not. allocated(outfield)) call error('You must specify field with -of option')
  fidx = 1 ! field number
  !do fidx=1,...
    pi = 1
    mtdim = 0
    do i=1, size(pmh%pc,1)
      if(allocated(outfield)) then
        if(size(outfield,1) /= get_piece_num_fields(pmh%pc(i))) &
          call error('Number of field file names and fields must agree')
      endif
      ! Point data
      if(allocated(pmh%pc(i)%fi)) then
        do j=1, size(pmh%pc(i)%fi,1)
          if(.true.) then ! Name control (?). Not used in mff
            filename = trim(path)//trim(outfield(fidx))
            if(.not. allocated(pmh%pc(i)%fi(j)%val)) &
               &call error("save_mff/ Point field "//trim(pmh%pc(i)%fi(j)%name)//": not allocated")
            call info('Writing node field "'//trim(adjustl(pmh%pc(i)%fi(j)%name))//'" to: '//trim(adjustl(filename)))
            if(present(param)) then
              do k=1, size(pmh%pc(i)%fi(j)%param,1)
                if((pmh%pc(i)%fi(j)%param(k)-param)<pmh%ztol) pi = k
              enddo
            endif
            iu = get_unit()
            open (unit=iu, file=trim(filename), form='formatted', position='rewind', iostat=ios)
            if (ios /= 0) call error('save/open, #'//trim(string(ios)))
            !write(unit=iu, fmt=*, iostat = ios) size(pmh%pc(i)%z,2)*size(pmh%pc(i)%fi(j)%val,1)
            call feed(iu, string(size(pmh%pc(i)%z,2)*size(pmh%pc(i)%fi(j)%val,1)))
            !if (ios /= 0) call error('save_mff/header, #'//trim(string(ios)))
            do k=1,size(pmh%pc(i)%fi(j)%val,2)
              !do l=1,size(pmh%pc(i)%fi(j)%val,1)
              !write(unit=iu, fmt=*, iostat = ios) &
              !   & (pmh%pc(i)%fi(j)%val(l,k,pi), l=1, size(pmh%pc(i)%fi(j)%val,1) )
              do l=1, size(pmh%pc(i)%fi(j)%val,1)
                call feed(iu, string(pmh%pc(i)%fi(j)%val(l,k,pi)))
              end do
              !if (ios /= 0) call error('save_mff/header, #'//trim(string(ios)))
              !enddo
            enddo
            call empty(iu)
            close(iu)
            fidx = fidx + 1
          endif
        enddo
      endif
      ! Cell data
      do j=1, size(pmh%pc(i)%el,1)
        mtdim = max(FEDB(pmh%pc(i)%el(j)%type)%tdim,mtdim)
      enddo
      do j=1, size(pmh%pc(i)%el,1)
        if (mtdim == FEDB(pmh%pc(i)%el(j)%type)%tdim .and. allocated(pmh%pc(i)%el(j)%fi)) then
          do k=1,size(pmh%pc(i)%el(j)%fi,1)
            if(.true.) then ! Name control(?). Not used in mff
              filename = trim(path)//trim(outfield(fidx))
              if(.not. allocated(pmh%pc(i)%el(j)%fi(k)%val)) &
                & call error("save_mff/ Cell field "//trim(pmh%pc(i)%el(j)%fi(k)%name)//": not allocated")
              call info('Writing cell field "'//trim(adjustl(pmh%pc(i)%el(j)%fi(k)%name))//'" to: '//trim(adjustl(filename)))
              if(present(param)) then
                do l=1, size(pmh%pc(i)%el(j)%fi(l)%param,1)
                  if((pmh%pc(i)%el(j)%fi(k)%param(l)-param)<pmh%ztol) pi = k
                enddo
              endif
              ! write pmh%pc(i)%el(j)%nel
              iu = get_unit()
              open (unit=iu, file=trim(filename), form='formatted', position='rewind', iostat=ios)
              if (ios /= 0) call error('save_mff/open, #'//trim(string(ios)))
              call feed(iu, string(pmh%pc(i)%el(j)%nel*size(pmh%pc(i)%el(j)%fi(k)%val,1)))
              !write(unit=iu, fmt=*, iostat = ios) pmh%pc(i)%el(j)%nel*size(pmh%pc(i)%el(j)%fi(k)%val,1)
              !if (ios /= 0) call error('save_mff/header, #'//trim(string(ios)))
              do l=1,size(pmh%pc(i)%el(j)%fi(k)%val,2)
                !write(unit=iu, fmt=*, iostat = ios) &
                !& (pmh%pc(i)%el(j)%fi(k)%val(m,l,pi), m=1, size(pmh%pc(i)%el(j)%fi(k)%val,1) )
                !if (ios /= 0) call error('save_mff/header, #'//trim(string(ios)))
                do m = 1, size(pmh%pc(i)%el(j)%fi(k)%val,1)
                  call feed(iu, string(pmh%pc(i)%el(j)%fi(k)%val(m,l,pi)))
                enddo
              enddo
              call empty(iu)
              close(iu)
              fidx = fidx + 1
            endif
          enddo
        endif
      enddo
    enddo
  !enddo
  print '(a)', 'Done!'
end subroutine

subroutine fix_filename(filename)
  character(len=*), intent(inout) :: filename !file names
  character(1), dimension(9)      ::  chars = ['<','>',':','"','/','\','|','?','*']
  integer                         :: i

  do i=1, size(chars,1)
    call replace(filename,chars(i),'_')
  enddo
end subroutine

end module
