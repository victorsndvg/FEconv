module module_ip_fcnv
!-----------------------------------------------------------------------
! Module to manage IP (Interpolation points) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Víctor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/05/2014
!
! PUBLIC PROCEDURES:
!   load_ip: loads a field from a IP format file
!   save_ip: saves a field in a IP format file
!-----------------------------------------------------------------------
use basicmod, only: real64, get_unit, string, replace, error
use module_utils_msh_fcnv
use module_pmh_fcnv

implicit none!


contains

subroutine load_ip(pmh, filenames, infieldnames, outfieldnames, param)
  type tempfield
    real(real64), allocatable               :: val(:,:)
  endtype
  character(len=*), allocatable, intent(in) :: filenames(:) !fields file names
  type(pmh_mesh),             intent(inout) :: pmh
  character(*), allocatable,  intent(inout) :: infieldnames(:) !Inpùt field names
  character(*), allocatable,  intent(inout) :: outfieldnames(:) !Output field names
  real(real64), optional,        intent(in) :: param
  character(len=maxpath)                    :: filename, aux
  integer                                   :: iu, ios, i, j, k, idx, maxtdimtotnel
  integer                                   :: ncomp, maxtdim, compcount
  integer                                   :: version, n_points, n_comps, n_fields, ncells, dim
  integer, allocatable                      :: fieldcomp(:,:) ! [(field number, number of component) x every component]
  integer, allocatable                      :: compsperfield(:) ! [(field number, number of component) x every component]
  character(len=maxpath), allocatable       :: compnames(:), fnames(:) !comp names
  type(field), allocatable                  :: tempfields(:)
  logical, allocatable                      :: found_comp(:)
  type(tempfield), allocatable              :: tfields(:)



  if(.not. allocated(filenames)) call error('Filename must be specified to load IP files.')

  do j=1, size(filenames,1)

    if(size(pmh%pc,1) == 1) then
      filename = trim(filenames(j))

      iu = get_unit()

      call info('Reading fields from: '//trim(adjustl(filenames(j))))
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

      if(allocated(found_comp)) deallocate(found_comp)
      allocate(found_comp(n_comps))

      if(allocated(fnames)) deallocate(fnames)
      allocate(fnames(n_comps))

      if(allocated(fieldcomp)) deallocate(fieldcomp)
      allocate(fieldcomp(2,n_comps))

      if(allocated(compsperfield)) deallocate(compsperfield)
      allocate(compsperfield(n_comps))

      ! Search the max topological dimension
      maxtdim = 0
      do i=1, size(pmh%pc(1)%el,1)
        maxtdim = max(FEDB(pmh%pc(1)%el(i)%type)%tdim,maxtdim)
      enddo

      ! Count the total number of elements with max topological dimension
      maxtdimtotnel = 0
      do i=1, size(pmh%pc(1)%el,1)
        if(FEDB(pmh%pc(1)%el(i)%type)%tdim == maxtdim) &
          maxtdimtotnel = maxtdimtotnel + pmh%pc(1)%el(i)%nel
      enddo

!      ! Read fieldnames and searchs patterns to build vector fields
!      n_fields = 0
!      do i=1,n_comps
!        is_vector_comp = .false.
!        read(unit=iu, fmt=*, iostat=ios) compnames(i)
!        if (ios /= 0) call error('load_ip/fieldnames, #'//trim(string(ios)))
!        comp = [index(compnames(i),'x-'), index(compnames(i),'y-'), index(compnames(i),'z-')]
!        ncomp = maxloc(comp,1)
!        if(maxval(comp)==1 .and. (ncomp>=1 .or. ncomp<=3)) then ! x, y and z
!          do k=1, i-1
!            if(index(compnames(k),compnames(i)(3:))>0) then
!              fieldcomp(1,i) = fieldcomp(1,k) ! Field number
!              fieldcomp(2,i) = ncomp ! Component number
!              compsperfield(fieldcomp(1,k)) = max(fieldcomp(2,k),ncomp)
!              fnames(n_fields) = trim(compnames(i)(3:))
!              is_vector_comp = .true.
!              exit
!            endif
!          enddo
!          if(.not. is_vector_comp) then
!            n_fields = n_fields + 1
!            fieldcomp(1,i) = n_fields ! Field number
!            fieldcomp(2,i) = ncomp ! Component number
!            compsperfield(n_fields) = 1
!            fnames(n_fields) = trim(compnames(i))
!          endif
!        else
!          n_fields = n_fields + 1
!          fieldcomp(1,i) = n_fields ! Field number
!          fieldcomp(2,i) = ncomp ! Component number
!          compsperfield(n_fields) = 1
!          fnames(n_fields) = trim(compnames(i))
!        endif
!      enddo


      ! Read fieldnames and searchs patterns to build vector fields
      n_fields = 0
      compcount = 0
      ncomp = 0
      found_comp = .true.
      fieldcomp(:,:) = 0
      compsperfield(:) = 0
      do i=1,n_comps
        read(unit=iu, fmt=*, iostat=ios) compnames(i)
        if (ios /= 0) call error('load_ip/fieldnames, #'//trim(string(ios)))
        if(allocated(infieldnames)) then
          found_comp(i) = .false.
          do k=1, size(infieldnames,1)
            if(trim(compnames(i)) == trim(infieldnames(k))) found_comp(i)=.true.
          enddo
        endif

        if(found_comp(i)) then
          if(allocated(outfieldnames)) then
            if(compcount >= ncomp) then
              n_fields = n_fields + 1
              compcount = 1
            else
              compcount = compcount + 1
            endif

            if(size(outfieldnames,1)<n_fields) call error('Input components and output fieldnames not agree!')

            ncomp = get_field_components(outfieldnames(n_fields))
            fnames(n_fields) = get_field_name(outfieldnames(n_fields))
            compsperfield(n_fields) = ncomp
            fieldcomp(1,i) = n_fields
            fieldcomp(2,i) = compcount

          else
            n_fields = n_fields + 1
            fieldcomp(1,i) = n_fields
            fieldcomp(2,i) = 1
            compsperfield(n_fields) = 1
            fnames(n_fields) = trim(compnames(i))

          endif
        endif

      enddo

      ! Read coordinates for IP version 2 or 3
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

      ! Read field values
      do k=1,n_comps
        if(fieldcomp(1,k)<=0) cycle
        if(.not. allocated(tfields(fieldcomp(1,k))%val)) &
          & allocate(tfields(fieldcomp(1,k))%val(compsperfield(n_fields),n_points))
        i = 1
        do
          if(i > n_points) exit
          read(iu,fmt='(A)', iostat=ios) aux
          if (ios /= 0) call error('load_ip/values, #'//trim(string(ios)))
          call replace_char(aux, '(', ' ')
          call replace_char(aux, ')', ' ')
          if(is_blank_line(aux)) cycle
          if(fieldcomp(1,k)==0) then
            read(aux,fmt=*,iostat=ios) aux
          else
            read(aux,fmt=*,iostat=ios) tfields(fieldcomp(1,k))%val(fieldcomp(2,k),i)
          endif
          i = i + 1
        enddo
      enddo

      close(iu)


      ! Assign fields to PMH Structure
      do k=1,n_fields
        ncells = 0
        ! Field over nodes. 1,2 or 3 components allowed
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
          call info('  Assigning node field: '//trim(fnames(k)))
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
          pmh%pc(1)%fi(idx)%val(:,:,1) = tfields(k)%val(:compsperfield(k),:)
        else
          ! Field over cells. 1,2 or 3 components allowed
          do i=1, size(pmh%pc(1)%el,1)
            if(maxtdim == FEDB(pmh%pc(1)%el(i)%type)%tdim) then
              ! Max top dim elements can be in different element grops
              if(n_points == pmh%pc(1)%el(i)%nel .or. n_points == maxtdimtotnel) then
                if(.not. allocated(pmh%pc(1)%el(i)%fi)) then
                  allocate(pmh%pc(1)%el(i)%fi(1))
                else
                  if(allocated(tempfields)) deallocate(tempfields)
                  allocate(tempfields(size(pmh%pc(1)%el(i)%fi,1)+1))
                  tempfields(1:size(pmh%pc(1)%el(i)%fi,1)) = pmh%pc(1)%el(i)%fi(:)
                  call move_alloc(from=tempfields, to=pmh%pc(1)%el(i)%fi)
                endif
                idx = size(pmh%pc(1)%el(i)%fi,1)
                call info('  Assigning cell field: '//trim(fnames(k)))
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
                pmh%pc(1)%el(i)%fi(idx)%val(:,:,1) = tfields(k)%val(:compsperfield(k),ncells+1:ncells+pmh%pc(1)%el(i)%nel)
                ncells = ncells+pmh%pc(1)%el(i)%nel
                exit
              endif
            endif
          enddo
        endif
      enddo
    endif
  enddo
  print*, ''

  if(allocated(outfieldnames) .and. allocated(infieldnames)) then
    deallocate(infieldnames)
    allocate(infieldnames(size(outfieldnames,1)))
    do i=1, size(outfieldnames,1)
      outfieldnames(i) = get_field_name(outfieldnames(i))
      infieldnames(i) = outfieldnames(i)
    enddo
  endif

end subroutine


subroutine save_ip(pmh, filenames, infieldnames, outfieldnames, nparam)
  type tempfield
    real(real64), allocatable               :: val(:,:)
  endtype
  character(len=*), allocatable, intent(in) :: filenames(:) !fields file names
  type(pmh_mesh),             intent(inout) :: pmh
  character(*), allocatable,  intent(inout) :: infieldnames(:) !Output field names
  character(*), allocatable,  intent(inout) :: outfieldnames(:) !Output field names
  integer, optional,             intent(in) :: nparam
  character(len=maxpath)                    :: filename, str,inichar
  integer :: maxdim, maxtopdim, n_points, n_comps, counter, lnv
  integer :: i, j, k, l, m, n, o, p, iu, ios, np
  integer, allocatable :: piece2save(:)
  logical :: node_fields, cell_fields, all_P1, check
  real(real64), allocatable :: znod(:,:), cellcoords(:,:)


  np = 1
  if(present(nparam)) np = nparam


  !check piece(s) to be saved
  if (is_arg('-p')) then !save a single piece, indicated after -p
    str = get_post_arg('-p')
    call alloc(piece2save, 1)
    piece2save(1) = int(str)
  else !save all pieces
    call alloc(piece2save, size(pmh%pc,1))
    piece2save = [(i, i=1, size(pmh%pc,1))]
  end if

  if(size(piece2save,1)>1) call error('save_ip/ Only one piece allowed')

  maxdim = 0
  n_points = 0
  node_fields = .false.
  cell_fields = .false.
  do i=1,size(piece2save,1)
    maxdim = max(maxdim, pmh%pc(piece2save(i))%dim)
    maxtopdim = get_piece_max_top_dim(pmh%pc(piece2save(i)))
    if(get_piece_num_fields(pmh%pc(piece2save(i))) == 0) call error('Piece '//trim(string(i))//' with no fields')
    if(allocated(infieldnames)) then
      n_comps = 0
      if(allocated(pmh%pc(piece2save(i))%fi)) then
        check = .false.
        do j=1, size(pmh%pc(piece2save(i))%fi,1)
          do k=1, size(infieldnames,1)
            if(trim(pmh%pc(piece2save(i))%fi(j)%name)==trim(infieldnames(k))) then
              check = .true.
              node_fields = .true.
              n_comps = n_comps + size(pmh%pc(piece2save(i))%fi(j)%val,1)
            endif
          enddo
        enddo
        if(check) n_points = n_points + pmh%pc(piece2save(i))%nnod
      endif
      counter = 1
      do j=1, size(pmh%pc(piece2save(i))%el,1)
        check = .false.
        if(maxtopdim /= FEDB(pmh%pc(piece2save(i))%el(j)%type)%tdim) cycle
        if(allocated(pmh%pc(piece2save(i))%el(j)%fi)) then
          n_comps = 0
          do k=1, size(pmh%pc(piece2save(i))%el(j)%fi,1)
            do l=1, size(infieldnames,1)
              if(trim(pmh%pc(piece2save(i))%el(j)%fi(k)%name)==trim(infieldnames(l))) then
                check = .true.
                cell_fields = .true.
                n_comps = n_comps + size(pmh%pc(piece2save(i))%el(j)%fi(k)%val,1)
              endif
            enddo
          enddo
          if(check) then
            ! Cell centroid coordinates
            n_points = n_points + pmh%pc(piece2save(i))%el(j)%nel
            lnv = FEDB(pmh%pc(piece2save(i))%el(j)%type)%lnv
            do k=1,pmh%pc(piece2save(i))%el(j)%nel
              do l=1,lnv
                call add(2, cellcoords, &
                       & pmh%pc(piece2save(i))%z(1:maxdim,pmh%pc(piece2save(i))%el(j)%mm(l,k)), &
                       &  counter, fit=[.true.,.false.])
              enddo
              cellcoords(1:maxdim,counter) = cellcoords(1:maxdim,counter)/lnv
              counter = counter + 1
            enddo
            call reduce(cellcoords, maxdim, n_points)
          endif
        endif
      enddo

    else

      if(get_piece_num_fields(pmh%pc(piece2save(i)),'node') == 0) then
        cell_fields = .true.
        counter = 1
        do j=1, size(pmh%pc(piece2save(i))%el,1)
          if(maxtopdim /= FEDB(pmh%pc(piece2save(i))%el(j)%type)%tdim) cycle
          n_points = n_points + pmh%pc(piece2save(i))%el(j)%nel
          n_comps = 0
          if(allocated(pmh%pc(piece2save(i))%el(j)%fi)) then
            do k=1, size(pmh%pc(piece2save(i))%el(j)%fi,1)
              n_comps = n_comps + size(pmh%pc(piece2save(i))%el(j)%fi(k)%val,1)
            enddo
            ! Cell centroid coordinates
            lnv = FEDB(pmh%pc(piece2save(i))%el(j)%type)%lnv
            do k=1,pmh%pc(piece2save(i))%el(j)%nel
              do l=1,lnv
                call add(2, cellcoords, &
                       & pmh%pc(piece2save(i))%z(1:maxdim,pmh%pc(piece2save(i))%el(j)%mm(l,k)), &
                       &  counter, fit=[.true.,.false.])
              enddo
              cellcoords(1:maxdim,counter) = cellcoords(1:maxdim,counter)/lnv
              counter = counter + 1
            enddo
          endif
        enddo
        call reduce(cellcoords, maxdim, n_points)
      elseif(get_piece_num_fields(pmh%pc(piece2save(i)),'cell') == 0) then
        node_fields = .true.
        n_points = n_points + pmh%pc(piece2save(i))%nnod
        n_comps = 0
        if(allocated(pmh%pc(piece2save(i))%fi)) then
          do j=1, size(pmh%pc(piece2save(i))%fi,1)
            n_comps = n_comps + size(pmh%pc(piece2save(i))%fi(j)%val,1)
          enddo
        endif
      else
        call error('save_ip/Cannot save node and cell fields simultaneously')
      endif
    endif
  enddo

  if(.not. allocated(filenames)) call error('save_ip/Output filename must be specified.')
  if(size(filenames,1)>1) call error('save_ip/Multiple files not allowed.')

  if(node_fields .and. cell_fields) call error('save_ip/Cannot save node and cell fields simultaneously')

  if(allocated(outfieldnames)) then
    if(n_comps /= size(outfieldnames,1)) call error('save_ip/Number of field names and component not agree.')
  endif

  do l=1, size(filenames,1)
    filename = trim(filenames(l))
    iu = get_unit()

    call info('Saving fields to: '//trim(adjustl(filenames(l))))
    !open file
    open (unit=iu, file=filename, form='formatted', status='replace', position='rewind', action="write",iostat=ios)
    if (ios /= 0) call error('save_ip/open, #'//trim(string(ios)))

    write(unit=iu, fmt= '(A)', iostat=ios) trim(string(3))
    if (ios /= 0) call error('save_ip/version, #'//trim(string(ios)))

    write(unit=iu, fmt= '(A)', iostat=ios) trim(string(maxdim))
    if (ios /= 0) call error('save_ip/dim, #'//trim(string(ios)))

    write(unit=iu, fmt= '(A)', iostat=ios) trim(string(n_points))
    if (ios /= 0) call error('save_ip/n_points, #'//trim(string(ios)))

    write(unit=iu, fmt='(A)', iostat=ios) trim(string(n_comps))
    if (ios /= 0) call error('save_ip/n_comps, #'//trim(string(ios)))

    if(allocated(outfieldnames)) then
      do i=1, size(outfieldnames,1)
        write(unit=iu, fmt='(A)', iostat=ios) trim(outfieldnames(i))
        if (ios /= 0) call error('save_ip/Component names, #'//trim(string(ios)))
      enddo
    else
      do i=1,n_comps
        write(unit=iu, fmt='(A)', iostat=ios) 'uds-'//trim(string(i-1))
        if (ios /= 0) call error('save_ip/Component names (uds), #'//trim(string(ios)))
      enddo
    endif

    do j=1,maxdim
      inichar = '('

      do i=1, size(piece2save,1)
        if(node_fields) then
          call build_node_coordinates(pmh%pc(piece2save(i)), i, all_P1, znod)
          if(all_P1) then
            if(size(pmh%pc(piece2save(i))%z,2)<n_points) call error('save_ip/Wrong number of point coordinates')
            do k=1, size(pmh%pc(piece2save(i))%z,2)
              write(unit=iu, fmt='(A)', iostat=ios) trim(inichar)//trim(string(pmh%pc(piece2save(i))%z(j,k)))
              if (ios /= 0) call error('save_ip/Z Coordinates, #'//trim(string(ios)))
              if(k==size(pmh%pc(piece2save(i))%z,2)) write(unit=iu, fmt='(A)', iostat=ios) ')'
              inichar = ''
            enddo
          else
            if(size(znod,2)<n_points) call error('save_ip/Wrong number of point coordinates')
            do k=1, size(znod,2)
              write(unit=iu, fmt='(A)', iostat=ios) trim(inichar)//trim(string(pmh%pc(piece2save(i))%z(j,k)))
              if (ios /= 0) call error('save_ip/Znod Coordinates, #'//trim(string(ios)))
              if(k==size(znod,2)) write(unit=iu, fmt='(A)', iostat=ios) ')'
              inichar = ''
            enddo
          endif
          inichar = ''
        elseif(cell_fields) then
          if(allocated(cellcoords)) then
            if(size(cellcoords,2)<n_points) call error('save_ip/Wrong number of cell coordinates')
            do k=1, size(cellcoords,2)
              write(unit=iu, fmt='(A)', iostat=ios) trim(inichar)//trim(string(cellcoords(j,k)))
              if (ios /= 0) call error('save_ip/Cell centroid coordinates, #'//trim(string(ios)))
              if(k==size(cellcoords,2)) write(unit=iu, fmt='(A)', iostat=ios) ')'
              inichar = ''
            enddo
          else
            call error('save_ip/Cell centroid coordinates not allocated')
          endif
        endif
      enddo
    enddo

    do i=1, size(piece2save,1)
      if(node_fields) then
        do j=1, size(pmh%pc(piece2save(i))%fi,1)
          do k=1, size(pmh%pc(piece2save(i))%fi(j)%val,1)
            inichar = '('
            do m=1, size(pmh%pc(piece2save(i))%fi(j)%val,2)
              write(unit=iu, fmt='(A)', iostat=ios) trim(inichar)//trim(string(pmh%pc(piece2save(i))%fi(j)%val(k,m,np)))
              if (ios /= 0) call error('save_ip/Node values, #'//trim(string(ios)))
              if(m==size(pmh%pc(piece2save(i))%fi(j)%val,2) ) write(unit=iu, fmt='(A)', iostat=ios) ')'
              inichar = ''
            enddo
          enddo
        enddo
      elseif(cell_fields) then
        do k=1, size(pmh%pc(piece2save(i))%el,1)
          check = .false.
          if(maxtopdim /= FEDB(pmh%pc(piece2save(i))%el(k)%type)%tdim) cycle
          if(allocated(pmh%pc(piece2save(i))%el(k)%fi)) then
            n_comps = 0
            do p=1, size(pmh%pc(piece2save(i))%el(k)%fi,1)
              inichar = '('
              if(allocated(infieldnames)) then
                do m=1, size(infieldnames,1)
                  if(trim(pmh%pc(piece2save(i))%el(k)%fi(p)%name)==trim(infieldnames(m))) then
                    do n=1,size(pmh%pc(piece2save(i))%el(k)%fi(p)%val,1)
                      do o=1,size(pmh%pc(piece2save(i))%el(k)%fi(p)%val,2)
                        write(unit=iu, fmt='(A)', iostat=ios) &
                          & trim(inichar)//trim(string(pmh%pc(piece2save(i))%el(k)%fi(p)%val(n,o,np)))
                        if (ios /= 0) call error('save_ip/Node values, #'//trim(string(ios)))
                        if(o==size(pmh%pc(piece2save(i))%el(k)%fi(p)%val,2)) write(unit=iu, fmt='(A)', iostat=ios) ')'
                        inichar = ''
                      enddo
                    enddo
                  endif
                enddo
              else
                if(allocated(pmh%pc(piece2save(i))%el(k)%fi(p)%val)) then
                  do n=1,size(pmh%pc(piece2save(i))%el(k)%fi(p)%val,1)
                    do o=1,size(pmh%pc(piece2save(i))%el(k)%fi(p)%val,2)
                      write(unit=iu, fmt='(A)', iostat=ios) &
                        & trim(inichar)//trim(string(pmh%pc(piece2save(i))%el(k)%fi(p)%val(n,o,np)))
                      if (ios /= 0) call error('save_ip/Node values, #'//trim(string(ios)))
                      if(o==size(pmh%pc(piece2save(i))%el(k)%fi(p)%val,2)) write(unit=iu, fmt='(A)', iostat=ios) ')'
                      inichar = ''
                    enddo
                  enddo
                endif
              endif
            enddo
          endif
        enddo
      endif
    enddo

  enddo
  print*, ''


end subroutine



subroutine fix_filename(filename)
  character(len=*), intent(inout) :: filename !file names
  character(1), dimension(9)      ::  chars = ['<','>',':','"','/','\','|','?','*']
  integer                         :: i

  do i=1, size(chars,1)
    call replace(filename,chars(i),'_')
  enddo

end subroutine

function get_field_components(fieldname) result(res)
  character(len=*), intent(in) :: fieldname
  integer                      :: res, i1,i2, ios

  res = 0
  if(len_trim(fieldname) == 0) then
    res = 0
  else
    i1 = index(fieldname,'(')
    i2 = index(fieldname,')')
    if(len_trim(fieldname) == 0) then
      call error('Field without name!')
    elseif((i1==0 .and. i2/=0) .or. (i1/=0 .and. i2==0)) then
      call error('Parentheses can appear only in the field names to indicate the number of components')
    elseif(i1 == 0) then
      res = 1
    elseif(i2>i1 .and. i2-i1>=2) then
      read(fieldname(i1+1:i2-1),fmt=*,iostat=ios) res
    else
      call error('Cannot parse number of componentes in field "'//trim(fieldname)//'"')
    endif
  endif

end function

function get_field_name(fieldname) result(res)
  character(len=*), intent(in) :: fieldname
  character(len=maxpath)       :: res
  integer                      :: i1,i2

  res = trim(fieldname)
  if(len_trim(fieldname) == 0) then
    res = trim(fieldname)
  else
    i1 = index(fieldname,'(')
    i2 = index(fieldname,')')
    if(len_trim(fieldname) == 0) then
      call error('Field without name!')
    elseif((i1==0 .and. i2/=0) .or. (i1/=0 .and. i2==0)) then
      call error('Parentheses can appear only in the field names to indicate the number of components')
    elseif(i1 == 0) then
      res = trim(fieldname)
    elseif(i2>i1 .and. i2-i1>=2) then
      res = trim(fieldname(1:i1-1))
    else
      call error('Cannot parse number of componentes in field "'//trim(fieldname)//'"')
    endif
  endif

end function

end module
