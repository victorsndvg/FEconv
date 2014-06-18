module module_dex
!-----------------------------------------------------------------------
! Module to manage DEX (Flux field) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Víctor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/05/2014
!
! PUBLIC PROCEDURES:
!   load_dex: loads a mesh from a DEX format file
!   save_dex: saves a mesh in a DEX format file
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

subroutine load_dex(pmh, filenames, fieldnames, param)
  character(len=*), allocatable, intent(in) :: filenames(:) !fields file names
  type(pmh_mesh),             intent(inout) :: pmh
  character(*), allocatable,     intent(in) :: fieldnames(:) !field names
  real(real64), optional,        intent(in) :: param 
  character(len=maxpath)                    :: filename, fieldname, aux
  integer                                   :: nb_real, nb_comp, nb_point
  integer                                   :: iu, ios, i, j, idx
  real(real64), allocatable                 :: coords(:,:), vals(:,:)
  type(field), allocatable                  :: tempfields(:)
!  integer                                   :: ncomp, totcomp, maxtdim



  if(size(filenames,1) /= size(fieldnames,1)) &
    call error('load_dex/ Filenames and fieldnames dimension must agree')

  do j=1, size(filenames,1)
    filename = trim(filenames(j))
    fieldname = trim(fieldnames(j))
    iu = get_unit()
    !open file

    open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
    if (ios /= 0) call error('load/open, #'//trim(string(ios)))
  
    !try read field name
    read(unit=iu, fmt=*, iostat=ios) aux,aux,aux,aux,aux,aux,aux
    if (ios /= 0) call error('load/open, #'//trim(string(ios)))
    !try read number of real, number of components and number of points
    read(unit=iu, fmt=*, iostat=ios) aux,aux,nb_real,aux,aux,nb_comp,aux,aux,nb_point
    if (ios /= 0) call error('load/open, #'//trim(string(ios)))  

    if(size(pmh%pc,1) == 1) then
      if(pmh%pc(1)%nnod /= nb_point) then
        call info('Number of values in field must agree with number of nodes. Skipped!')
      else
        if(allocated(coords)) deallocate(coords)
        if(allocated(vals)) deallocate(vals)
        allocate(coords(pmh%pc(1)%dim,nb_point))
        allocate(vals(nb_comp,nb_point))
        ! Read coords and values
        do i=1,nb_point
          read(unit=iu, fmt=*, iostat=ios) coords(:,i),vals(:,i)
          if (ios /= 0) call error('load/open, #'//trim(string(ios)))  
        enddo
        if(.not. allocated(pmh%pc(1)%fi)) then 
          allocate(pmh%pc(1)%fi(1))
        else
          if(allocated(tempfields)) deallocate(tempfields)
          allocate(tempfields(size(pmh%pc(1)%fi,1)+1))
          tempfields(1:size(pmh%pc(1)%fi,1)) = pmh%pc(1)%fi(:)
          call move_alloc(from=tempfields, to=pmh%pc(1)%fi)
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
        allocate(pmh%pc(1)%fi(idx)%val(nb_comp, nb_point,1))
        pmh%pc(1)%fi(idx)%val(:,:,1) = vals(:,:)

      endif
    endif

    close(iu)

  enddo
  print*, ''

end subroutine


subroutine save_dex(pmh, infield, outfield, path, param)
  type(pmh_mesh),            intent(inout) :: pmh      !PMH mesh
  character(*), allocatable, intent(in) :: infield(:)  ! In field names
  character(*), allocatable, intent(in) :: outfield(:) ! Out field names
  character(*),              intent(in) :: path !file names
  real(real64), optional,    intent(in) :: param 
  real(real64),allocatable              :: znod(:,:)
  character(len=maxpath)                :: filename !file names
  integer                               :: i,j,k,l,pi,mtdim
  integer                               :: iu, ios, fidx
  logical                               :: all_f, all_P1

  if(size(infield,1) /= size(outfield,1)) &
    call error('load_dex/ Filenames and fieldnames dimension must agree')

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
        if(allocated(znod)) deallocate(znod)
        call build_node_coordinates(pmh%pc(i), i, all_P1, znod)
        do j=1, size(pmh%pc(i)%fi,1)
          if(trim(infield(fidx)) == trim(pmh%pc(i)%fi(j)%name) .or. all_f) then
            if(.not. allocated(pmh%pc(i)%fi(j)%val)) &
               &call error("save_dex/ Point field "//trim(infield(fidx))//": not allocated")
            call fix_filename(pmh%pc(i)%fi(j)%name)
            if(all_f .and. trim(infield(1)) == '*') &
              & filename = trim(path)//trim(outfield(1))//'__'//trim(pmh%pc(i)%fi(j)%name)//'.dex'
            call info('Writing node field to: '//trim(adjustl(filename)))
            if(present(param)) then
              do k=1, size(pmh%pc(i)%fi(j)%param,1)
                if((pmh%pc(i)%fi(j)%param(k)-param)<pmh%ztol) pi = k
              enddo
            endif
            if(pmh%pc(i)%nnod /= size(znod,2) .and. pmh%pc(i)%nnod /= size(pmh%pc(i)%fi(j)%val,2)) then
              call info('  Wrong number of values. Skipped!')
              cycle
            endif
            iu = get_unit() 
            open (unit=iu, file=trim(filename), form='formatted', position='rewind', iostat=ios)
            if (ios /= 0) call error('save_dex/header, #'//trim(string(ios)))

            ! Header
            write(unit=iu, fmt=*, iostat = ios) &
              & '# NAME = PIECE'//trim(string(i))// ' FORMULA = '//trim(adjustl(pmh%pc(i)%fi(j)%name))
            if (ios /= 0) call error('save_dex/header, #'//trim(string(ios)))
            write(unit=iu, fmt=*, iostat = ios) &
              & 'NB_REAL = '//trim(string(1))// ' NB_COMP = '//trim(string(size(pmh%pc(i)%fi(j)%val,1)))//&
              & ' NB_POINT = '//trim(string(size(pmh%pc(i)%fi(j)%val,2)))//' #'
            if (ios /= 0) call error('save_dex/header, #'//trim(string(ios)))

            if(all_P1) then
              do k=1,pmh%pc(i)%nver
                write(unit=iu, fmt=*, iostat = ios) &
                  & (pmh%pc(i)%z(l,k), l=1, size(pmh%pc(i)%z,1) ), &
                  & (pmh%pc(i)%fi(j)%val(l,k,pi), l=1, size(pmh%pc(i)%fi(j)%val,1) )
                if (ios /= 0) call error('save_dex/vertex, #'//trim(string(ios)))
              enddo
            else
              do k=1,size(znod,2)
                write(unit=iu, fmt=*, iostat = ios) &
                  & (znod(l,k), l=1, size(znod,1) ), &
                  & (pmh%pc(i)%fi(j)%val(l,k,pi), l=1, size(pmh%pc(i)%fi(j)%val,1) )
                if (ios /= 0) call error('save_dex/nodes, #'//trim(string(ios)))
              enddo
            endif
            close(iu)
          endif
        enddo
      endif

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

