module module_mff
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
use module_compiler_dependant, only: real64
use module_files, only: get_unit
use module_convers, only: string, replace
use module_report, only:error
use module_pmh

implicit none!

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
  character(*), allocatable,     intent(in) :: fieldnames(:) !field names
  real(real64), optional,        intent(in) :: param 
  character(len=maxpath)                    :: filename, fieldname
  integer                                   :: iu, ios, i, j, k, idx  
  integer                                   :: ncomp, totcomp, maxtdim
  real(real64), allocatable                 :: fielddata(:)
  type(field), allocatable                  :: tempfields(:)


  if(size(filenames,1) /= size(fieldnames,1)) &
    call error('load_mff/ Filenames and fieldnames dimension must agree')

  do j=1, size(filenames,1)
    filename = trim(filenames(j))
    fieldname = trim(fieldnames(j))
    iu = get_unit()
    !open file
    open (unit=get_unit(), file=filename, form='formatted', status='old', position='rewind', iostat=ios)
    if (ios /= 0) call error('load/open, #'//trim(string(ios)))
  
    !try read number of components
    read(unit=iu, fmt=*, iostat=ios) totcomp
  
    if(allocated(fielddata)) deallocate(fielddata)
    allocate(fielddata(totcomp))
  
    read(unit=iu, fmt=*, iostat=ios) (fielddata(i),  i=1,totcomp)
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
          deallocate(tempfields)
        endif
        idx = size(pmh%pc(1)%fi,1)
        call info('Reading node field from: '//trim(adjustl(filenames(j))))
        pmh%pc(1)%fi(idx)%name = trim(filename)
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
              call info('Reading cell field from: '//trim(adjustl(filenames(j))))
              pmh%pc(1)%el(i)%fi(idx)%name = trim(filename)
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
  print*, ''

end subroutine


subroutine save_mff(pmh, infield, outfield, path, param)
  type(pmh_mesh),            intent(inout) :: pmh      !PMH mesh
  character(*), allocatable, intent(in) :: infield(:)  ! In field names
  character(*), allocatable, intent(in) :: outfield(:) ! Out field names
  character(*),              intent(in) :: path !file names
  real(real64), optional,    intent(in) :: param 
  character(len=maxpath)                :: filename !file names
  integer                               :: i,j,k,l,m,pi,mtdim
  integer                               :: iu, ios, fidx
  logical                               :: all_f

  if(size(infield,1) /= size(outfield,1)) &
    call error('load_mff/ Filenames and fieldnames dimension must agree')

  all_f = .false.

  if(size(infield,1) == 1) all_f = (trim(infield(1)) == '*')
  if(size(infield,1) == 1 .and. size(outfield,1) == 1) all_f = .true.

  do fidx=1, size(infield,1)
    filename = trim(path)//trim(outfield(fidx))

    pi = 1
    mtdim = 0

    do i=1, size(pmh%pc,1)
      ! Point data
      if(allocated(pmh%pc(i)%fi)) then
        do j=1, size(pmh%pc(i)%fi,1)
          if(trim(infield(fidx)) == trim(pmh%pc(i)%fi(j)%name) .or. all_f) then
            if(.not. allocated(pmh%pc(i)%fi(j)%val)) &
               &call error("save_mff/ Point field "//trim(infield(fidx))//": not allocated")
            call fix_filename(pmh%pc(i)%fi(j)%name)
            if(all_f .and. trim(infield(1)) == '*') &
              & filename = trim(path)//trim(outfield(1))//'__'//trim(pmh%pc(i)%fi(j)%name)//'.mff'
            call info('Writing node field to: '//trim(adjustl(filename)))
            if(present(param)) then
              do k=1, size(pmh%pc(i)%fi(j)%param,1)
                if((pmh%pc(i)%fi(j)%param(k)-param)<pmh%ztol) pi = k
              enddo
            endif
            iu = get_unit() 
            open (unit=iu, file=trim(filename), form='formatted', position='rewind', iostat=ios)
            if (ios /= 0) call error('save/open, #'//trim(string(ios)))
            write(unit=iu, fmt=*, iostat = ios) size(pmh%pc(i)%z,2)*size(pmh%pc(i)%fi(j)%val,1)
            if (ios /= 0) call error('save_mff/header, #'//trim(string(ios)))
            do k=1,size(pmh%pc(i)%fi(j)%val,2)
  !            do l=1,size(pmh%pc(i)%fi(j)%val,1)
                write(unit=iu, fmt=*, iostat = ios) &
                  & (pmh%pc(i)%fi(j)%val(l,k,pi), l=1, size(pmh%pc(i)%fi(j)%val,1) )
                if (ios /= 0) call error('save_mff/header, #'//trim(string(ios)))
  !            enddo
            enddo
            close(iu)
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
            if(trim(infield(fidx)) == trim(pmh%pc(i)%el(j)%fi(k)%name)) then
              if(.not. allocated(pmh%pc(i)%el(j)%fi(k)%val)) &
                & call error("save_mff/ Cell field "//trim(infield(fidx))//": not allocated")
              call fix_filename(pmh%pc(i)%el(j)%fi(k)%name)
              if(all_f .and. trim(infield(1)) == '*') &
                & filename = trim(path)//trim(outfield(1))//'__'//trim(pmh%pc(i)%el(j)%fi(k)%name)//'.mff'
            call info('Writing cell field to: '//trim(adjustl(filename)))
              if(present(param)) then
                do l=1, size(pmh%pc(i)%el(j)%fi(l)%param,1)
                  if((pmh%pc(i)%el(j)%fi(k)%param(l)-param)<pmh%ztol) pi = k
                enddo
              endif
              ! write pmh%pc(i)%el(j)%nel
              do l=1,size(pmh%pc(i)%el(j)%fi(k)%val,2)
                do m=1,size(pmh%pc(i)%el(j)%fi(k)%val,1)
                  !write pmh%pc(i)%el(j)%fi(k)%val(l,k,pi)
                enddo
              enddo            
            endif
          enddo
        endif
      enddo
    enddo

  
  enddo
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

