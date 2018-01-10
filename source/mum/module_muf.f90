module module_muf_fcnv
!-----------------------------------------------------------------------
! Module to manage MUF (Modulef Unformatted Field) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Víctor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/05/2014
!
! PUBLIC PROCEDURES:
!   load_muf: loads a mesh from a MUM format file
!   save_muf: saves a mesh in a MUM format file
!-----------------------------------------------------------------------
use basicmod, only: real64, get_unit, string, replace, error
use module_pmh_fcnv
implicit none

!type field
!  character(maxpath)        :: name
!  real(real64), allocatable :: param(:)   !nshot
!  real(real64), allocatable :: val(:,:,:) !ncomp x nnod x nshot
!end type
!private :: field

contains

subroutine load_muf(pmh, filenames, fieldnames, param)
  character(len=*), allocatable, intent(in) :: filenames(:) !fields file names
  type(pmh_mesh),             intent(inout) :: pmh
  character(*), allocatable,     intent(in) :: fieldnames(:) !field names
  real(real64), optional,        intent(in) :: param
  character(len=maxpath)                    :: filename, fieldname
  integer                                   :: iu, ios, i, j, k, idx
  integer                                   :: ncomp, totcomp, maxtdim
  real(real64), allocatable                 :: fielddata(:)
  type(field), allocatable                  :: tempfields(:)


  if(allocated(fieldnames) .and. size(filenames,1) /= size(fieldnames,1)) &
    call error('load_muf/ Filenames and fieldnames dimension must agree')

  do j=1, size(filenames,1)
    filename = trim(filenames(j))
    fieldname = trim(filenames(j))

    if(allocated(fieldnames)) fieldname = trim(fieldnames(j))
    iu = get_unit()
    !open file
    open (unit=iu, file=filename, form='unformatted', status='old', position='rewind', iostat=ios)
    if (ios /= 0) call error('load/open, #'//trim(string(ios)))

    !try read number of components
    read(unit=iu, iostat=ios) totcomp

    if(allocated(fielddata)) deallocate(fielddata)
    allocate(fielddata(totcomp))

    read(unit=iu, iostat=ios) (fielddata(i),  i=1,totcomp)
    close(iu)

    maxtdim = 0

    if(size(pmh%pc,1) == 1) then
      ! Field over nodes. 1,2 or 3 components allowed
      if (mod(totcomp, size(pmh%pc(1)%z,2)) == 0 .and. (totcomp/size(pmh%pc(1)%z,2)<= 3)) then
        ncomp = totcomp/size(pmh%pc(1)%z,2)
        if(.not. allocated(pmh%pc(1)%fi)) then
          allocate(pmh%pc(1)%fi(1))
        else
          if(allocated(tempfields)) deallocate(tempfields)
          allocate(tempfields(size(pmh%pc(1)%fi,1)+1))
          tempfields(1:size(pmh%pc(1)%fi,1)) = pmh%pc(1)%fi(:)
          call move_alloc(from=tempfields, to=pmh%pc(1)%fi)
        endif
        idx = size(pmh%pc(1)%fi,1)
        call info('Reading node field "'//trim(adjustl(fieldname))//'" from: '//trim(adjustl(filename)))
        pmh%pc(1)%fi(idx)%name = trim(fieldname)
        if(allocated(pmh%pc(1)%fi(idx)%param)) deallocate(pmh%pc(1)%fi(idx)%param)
        allocate(pmh%pc(1)%fi(idx)%param(1))
        if(present(param)) then
          pmh%pc(1)%fi(idx)%param(1) = param
        else
          pmh%pc(1)%fi(idx)%param(1) = 0._real64
        endif
        if(allocated(pmh%pc(1)%fi(idx)%val)) deallocate(pmh%pc(1)%fi(idx)%val)
        allocate(pmh%pc(1)%fi(idx)%val(ncomp, totcomp/ncomp,1))
        do i=1, totcomp/ncomp
           pmh%pc(1)%fi(idx)%val(:,i,1) = fielddata((i-1)*ncomp+1:i*ncomp)
        enddo
      else
        do i=1, size(pmh%pc(1)%el,1)
          maxtdim = max(FEDB(pmh%pc(1)%el(i)%type)%tdim,maxtdim)
        enddo

        do i=1, size(pmh%pc(1)%el,1)
          if(maxtdim == FEDB(pmh%pc(1)%el(i)%type)%tdim) then
            if(mod(totcomp,pmh%pc(1)%el(i)%nel) == 0 .and. totcomp/pmh%pc(1)%el(i)%nel<= 3) then
              ncomp = totcomp/pmh%pc(1)%el(i)%nel
              if(.not. allocated(pmh%pc(1)%el(i)%fi)) then
                allocate(pmh%pc(1)%el(i)%fi(1))
              else
                if(allocated(tempfields)) deallocate(tempfields)
                allocate(tempfields(size(pmh%pc(1)%el(i)%fi,1)+1))
                tempfields(1:size(pmh%pc(1)%fi,1)) = pmh%pc(1)%el(i)%fi(:)
                call move_alloc(from=tempfields, to=pmh%pc(1)%el(i)%fi)
                deallocate(tempfields)
              endif
              idx = size(pmh%pc(1)%el(i)%fi,1)
              call info('Reading cell field "'//trim(adjustl(fieldname))//'" from: '//trim(adjustl(filename)))
              pmh%pc(1)%el(i)%fi(idx)%name = trim(fieldname)
              if(allocated(pmh%pc(1)%el(i)%fi(idx)%param)) deallocate(pmh%pc(1)%el(i)%fi(idx)%param)
              allocate(pmh%pc(1)%el(i)%fi(idx)%param(1))
              if(present(param)) then
                pmh%pc(1)%el(i)%fi(idx)%param(1) = param
              else
                pmh%pc(1)%el(i)%fi(idx)%param(1) = 0._real64
              endif
              if(allocated(pmh%pc(1)%el(i)%fi(idx)%val)) deallocate(pmh%pc(1)%el(i)%fi(idx)%val)
              allocate(pmh%pc(1)%el(i)%fi(idx)%val(ncomp, totcomp/ncomp,1))
              do k=1, pmh%pc(1)%el(i)%nel
                pmh%pc(1)%el(i)%fi(idx)%val(:,k,1) = fielddata((k-1)*ncomp+1:k*ncomp)
              enddo
              exit
            endif
          endif
        enddo
      endif
    endif
  enddo

print '(a)', 'Done!'

end subroutine


subroutine save_muf(pmh, outfield, path, param)
  type(pmh_mesh),            intent(inout) :: pmh      !PMH mesh
  character(*), allocatable, intent(in) :: outfield(:) ! Out field file names
  character(*),              intent(in) :: path        !path
  real(real64), optional,    intent(in) :: param
  character(len=maxpath)                :: filename    !file names
  integer                               :: i,j,k,l,m,pi,mtdim
  integer                               :: iu, ios, fidx



   if(.not. allocated(outfield)) call error('You must specify field with -of option')

  fidx = 1 ! field number
!  do fidx=1,


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
          if(.true.) then ! Name control (?). Not used in muf
            filename = trim(path)//trim(outfield(fidx))
            if(.not. allocated(pmh%pc(i)%fi(j)%val)) &
               &call error("save_muf/ Point field "//trim(pmh%pc(i)%fi(j)%name)//": not allocated")
            call info('Writing node field "'//trim(adjustl(pmh%pc(i)%fi(j)%name))//'" to: '//trim(adjustl(filename)))
            if(present(param)) then
              do k=1, size(pmh%pc(i)%fi(j)%param,1)
                if((pmh%pc(i)%fi(j)%param(k)-param)<pmh%ztol) pi = k
              enddo
            endif
            iu = get_unit()
            open (unit=iu, file=trim(filename), form='unformatted', position='rewind', iostat=ios)
            if (ios /= 0) call error('save/open, #'//trim(string(ios)))
            write(unit=iu, iostat = ios) size(pmh%pc(i)%z,2)*size(pmh%pc(i)%fi(j)%val,1)
            if (ios /= 0) call error('save_muf/header, #'//trim(string(ios)))
            do k=1,size(pmh%pc(i)%fi(j)%val,2)
  !            do l=1,size(pmh%pc(i)%fi(j)%val,1)
                write(unit=iu, iostat = ios) &
                  & (pmh%pc(i)%fi(j)%val(l,k,pi), l=1, size(pmh%pc(i)%fi(j)%val,1) )
                if (ios /= 0) call error('save_muf/header, #'//trim(string(ios)))
  !            enddo
            enddo
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
        if(mtdim == FEDB(pmh%pc(i)%el(j)%type)%tdim .and. allocated(pmh%pc(i)%el(j)%fi)) then
          do k=1,size(pmh%pc(i)%el(j)%fi,1)
            if(.true.) then ! Name control(?). Not used in muf
              filename = trim(path)//trim(outfield(fidx))
              if(.not. allocated(pmh%pc(i)%el(j)%fi(k)%val)) &
                & call error("save_muf/ Cell field "//trim(pmh%pc(i)%el(j)%fi(k)%name)//": not allocated")
              call info('Writing cell field "'//trim(adjustl(pmh%pc(i)%el(j)%fi(k)%name))//'" to: '//trim(adjustl(filename)))
              if(present(param)) then
                do l=1, size(pmh%pc(i)%el(j)%fi(l)%param,1)
                  if((pmh%pc(i)%el(j)%fi(k)%param(l)-param)<pmh%ztol) pi = k
                enddo
              endif
              ! write pmh%pc(i)%el(j)%nel
              iu = get_unit()
              open (unit=iu, file=trim(filename), form='unformatted', position='rewind', iostat=ios)
              if (ios /= 0) call error('save_muf/open, #'//trim(string(ios)))
              write(unit=iu, iostat = ios) pmh%pc(i)%el(j)%nel*size(pmh%pc(i)%el(j)%fi(k)%val,1)
              if (ios /= 0) call error('save_muf/header, #'//trim(string(ios)))
              do l=1,size(pmh%pc(i)%el(j)%fi(k)%val,2)
!                do m=1,size(pmh%pc(i)%el(j)%fi(k)%val,1)
                write(unit=iu, iostat = ios) &
                  & (pmh%pc(i)%el(j)%fi(k)%val(m,l,pi), m=1, size(pmh%pc(i)%el(j)%fi(k)%val,1) )
                if (ios /= 0) call error('save_muf/header, #'//trim(string(ios)))
!                enddo
              enddo
              close(iu)
              fidx = fidx + 1
            endif
          enddo
        endif
      enddo
    enddo
!  enddo

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
