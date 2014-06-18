module module_ip
!-----------------------------------------------------------------------
! Module to manage IP (Interpolation points) files
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
use module_utils_msh
use module_pmh

implicit none!


contains

subroutine load_ip(pmh, filenames, fieldnames, param)
  type tempfield
    real(real64), allocatable               :: val(:,:)
  endtype
  character(len=*), allocatable, intent(in) :: filenames(:) !fields file names
  type(pmh_mesh),             intent(inout) :: pmh
  character(*), allocatable,     intent(in) :: fieldnames(:) !field names
  real(real64), optional,        intent(in) :: param 
  character(len=maxpath)                    :: filename, fieldname, aux
  integer                                   :: iu, ios, i, j, k, idx  
  integer                                   :: ncomp, maxtdim, comp(3)
  integer                                   :: version, n_points, n_comps, n_fields, ncells, dim
  integer, allocatable                      :: fieldcomp(:,:) ! [(field number, number of component) x every component]
  integer, allocatable                      :: compsperfield(:) ! [(field number, number of component) x every component]
  character(len=maxpath), allocatable       :: compnames(:), fnames(:) !comp names
  type(field), allocatable                  :: tempfields(:)
  logical                                   :: is_vector_comp
  type(tempfield), allocatable              :: tfields(:)

  if(size(filenames,1) /= size(fieldnames,1)) &
    call error('load_ip/ Filenames and fieldnames dimension must agree')

  do j=1, size(filenames,1)

    if(size(pmh%pc,1) == 1) then
      filename = trim(filenames(j))
      fieldname = trim(fieldnames(j))
      iu = get_unit()
      !open file
      open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
      if (ios /= 0) call error('load_ip/open, #'//trim(string(ios)))
    
      read(unit=iu, fmt=*, iostat=ios) version
      if (ios /= 0) call error('load_ip/version, #'//trim(string(ios)))
  
      read(unit=iu, fmt=*, iostat=ios) dim
      if (ios /= 0) call error('load_ip/dim, #'//trim(string(ios)))
  
      read(unit=iu, fmt=*, iostat=ios) n_points
      if (ios /= 0) call error('load_ip/n_points, #'//trim(string(ios)))
  
      read(unit=iu, fmt=*, iostat=ios) n_comps
      if (ios /= 0) call error('load/n_comps, #'//trim(string(ios)))
  
      if(allocated(compnames)) deallocate(compnames)
      allocate(compnames(n_comps))
  
      if(allocated(fnames)) deallocate(fnames)
      allocate(fnames(n_comps))
  
      if(allocated(fieldcomp)) deallocate(fieldcomp)
      allocate(fieldcomp(2,n_comps))
  
      if(allocated(compsperfield)) deallocate(compsperfield)
      allocate(compsperfield(n_comps))


      n_fields = 0
      do i=1,n_comps
        is_vector_comp = .false.
        read(unit=iu, fmt=*, iostat=ios) compnames(i)
        if (ios /= 0) call error('load_ip/fieldnames, #'//trim(string(ios)))
        comp = [index(compnames(i),'x-'), index(compnames(i),'y-'), index(compnames(i),'z-')]
        if(maxval(comp)==1 .and. (ncomp>=1 .or. ncomp<=3)) then ! x, y and z
          ncomp = maxloc(comp,1)
          do k=1, i-1
            if(index(compnames(k),compnames(i)(3:))>0) then
              fieldcomp(1,i) = fieldcomp(1,k) ! Field number
              fieldcomp(2,i) = ncomp ! Component number
              compsperfield(fieldcomp(1,k)) = max(fieldcomp(2,k),ncomp)
              fnames(n_fields) = trim(compnames(i)(3:))
              is_vector_comp = .true.
              exit
            endif
          enddo
          if(.not. is_vector_comp) then
            n_fields = n_fields + 1
            fieldcomp(1,i) = n_fields ! Field number
            fieldcomp(2,i) = ncomp ! Component number
            compsperfield(n_fields) = 1
            fnames(n_fields) = trim(compnames(i))
          endif
        else
          n_fields = n_fields + 1
          fieldcomp(1,i) = n_fields ! Field number
          fieldcomp(2,i) = ncomp ! Component number
          compsperfield(n_fields) = 1
          fnames(n_fields) = trim(compnames(i))
        endif
      enddo  

      do k=1, dim
        i = 1
        do
          if(i > n_points) exit
          read(iu,fmt='(A)', iostat=ios) aux
          if (ios /= 0) call error('load_ip/coordinates, #'//trim(string(ios)))
          call replace_char(aux, '(', ' ') 
          call replace_char(aux, ')', ' ')
          if(is_blank_line(aux)) cycle
          i = i + 1
        enddo
      enddo

      if(allocated(tfields)) then
        do k=1, size(tfields)
          if(allocated(tfields(k)%val)) deallocate(tfields(k)%val)
        enddo
        deallocate(tfields)
      endif
      allocate(tfields(n_fields))

      do k=1,n_comps
        if(.not. allocated(tfields(k)%val)) allocate(tfields(k)%val(compsperfield(n_fields),n_points))
        i = 1
        do
          if(i > n_points) exit
          read(iu,fmt='(A)', iostat=ios) aux
          if (ios /= 0) call error('load_ip/values, #'//trim(string(ios)))
          call replace_char(aux, '(', ' ') 
          call replace_char(aux, ')', ' ')
          if(is_blank_line(aux)) cycle
          read(aux,fmt=*,iostat=ios) tfields(fieldcomp(1,k))%val(fieldcomp(2,k),i)
          i = i + 1
        enddo
      enddo

      close(iu)

      maxtdim = 0
      do i=1, size(pmh%pc(1)%el,1)
        maxtdim = max(FEDB(pmh%pc(1)%el(i)%type)%tdim,maxtdim)
      enddo
  
      do k=1,n_fields
        ncells = 0
        ! Field over cells. 1,2 or 3 components allowed
        if (n_points == pmh%pc(1)%nnod) then
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
          pmh%pc(1)%fi(idx)%name = trim(fnames(k))
          if(allocated(pmh%pc(1)%fi(idx)%param)) deallocate(pmh%pc(1)%fi(idx)%param)
          allocate(pmh%pc(1)%fi(idx)%param(1))      
          if(present(param)) then 
            pmh%pc(1)%fi(idx)%param(1) = param
          else
            pmh%pc(1)%fi(idx)%param(1) = 0._real64
          endif
          if(allocated(pmh%pc(1)%fi(idx)%val)) deallocate(pmh%pc(1)%fi(idx)%val)
          allocate(pmh%pc(1)%fi(idx)%val(compsperfield(k), n_points,1))
          pmh%pc(1)%fi(idx)%val = 0._real64
        else
        ! Field over cells. 1,2 or 3 components allowed
          do i=1, size(pmh%pc(1)%el,1)
            if(maxtdim == FEDB(pmh%pc(1)%el(i)%type)%tdim) then
! Max top dim elements can be in different element grops
!              if(n_points == pmh%pc(1)%el(i)%nel) then
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
                pmh%pc(1)%el(i)%fi(idx)%name = trim(fnames(k))
                if(allocated(pmh%pc(1)%el(i)%fi(idx)%param)) deallocate(pmh%pc(1)%el(i)%fi(idx)%param)
                allocate(pmh%pc(1)%el(i)%fi(idx)%param(1))      
                if(present(param)) then 
                  pmh%pc(1)%el(i)%fi(idx)%param(1) = param
                else
                  pmh%pc(1)%el(i)%fi(idx)%param(1) = 0._real64
                endif
                if(allocated(pmh%pc(1)%el(i)%fi(idx)%val)) deallocate(pmh%pc(1)%el(i)%fi(idx)%val)
                allocate(pmh%pc(1)%el(i)%fi(idx)%val(compsperfield(k), pmh%pc(1)%el(i)%nel,1))
                pmh%pc(1)%el(i)%fi(idx)%val = 0._real64
                pmh%pc(1)%el(i)%fi(idx)%val(:,:,1) = tfields(fieldcomp(1,k))%val(:,ncells+1:ncells+pmh%pc(1)%el(i)%nel)
                ncells = ncells+pmh%pc(1)%el(i)%nel
                exit
!              endif
            endif
          enddo
        endif
      enddo
    endif
  enddo
  print*, ''

  ! Memory deallocation
  if(allocated(compnames)) deallocate(compnames)
  if(allocated(fnames)) deallocate(fnames)
  if(allocated(fieldcomp)) deallocate(fieldcomp)
  if(allocated(compsperfield)) deallocate(compsperfield)
  if(allocated(tfields)) then
    do k=1, size(tfields)
      if(allocated(tfields(k)%val)) deallocate(tfields(k)%val)
    enddo
    deallocate(tfields)
  endif


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

