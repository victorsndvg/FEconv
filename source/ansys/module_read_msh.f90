module module_read_msh_fcnv

!-----------------------------------------------------------------------
! Module to manage MSH (Ansys) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/03/2014
!
! PUBLIC PROCEDURES:
!-----------------------------------------------------------------------

use basicmod
use module_pmh_fcnv
use module_utils_msh_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! read_msh_comment(iu, line): Read MSH comment section
!-----------------------------------------------------------------------
! iu: msh file unit number
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_comment(iu, line)
  character(len=*), intent(in) :: line
  integer,          intent(in) :: iu
  character(len=MAXPATH)       :: comm
  integer                      :: nopch,nclch, ios

  nopch = count_delimiter(line,'(')
  nclch = count_delimiter(line,')')
  do
    if(nopch <= nclch) exit
    read(iu,fmt='(A)', iostat=ios) comm
    if(ios /= 0) exit
    nopch = nopch+count_delimiter(comm,'(')
    nclch = nclch+count_delimiter(comm,')')
  enddo

end subroutine


!-----------------------------------------------------------------------
! read_msh_header(iu, line): Read MSH header section
!-----------------------------------------------------------------------
! iu: msh file unit number
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_header(iu, line)
  character(len=*), intent(in) :: line
  integer,          intent(in) :: iu

  call read_msh_comment(iu,line)

end subroutine


!-----------------------------------------------------------------------
! read_msh_section(iu, line): Read MSH generic section
!-----------------------------------------------------------------------
! iu: msh file unit number
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_section(iu, line)
  character(len=*), intent(in) :: line
  integer,          intent(in) :: iu

  call read_msh_comment(iu,line)

end subroutine


!-----------------------------------------------------------------------
! read_msh_dimensions(iu, line): Read MSH dimensions section
!-----------------------------------------------------------------------
! iu: msh file unit number
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_dimensions(pmh, line)
  character(len=*),       intent(in) :: line
  type(pmh_mesh),      intent(inout) :: pmh ! pmh_mesh
  character(len=MAXPATH)             :: auxline
  integer, allocatable               :: data(:)

  auxline(:) = line(:)
  call replace_char(auxline, '(', ' ')
  call replace_char(auxline, ')', ' ')
  if(allocated(data)) deallocate(data)
  !allocate(data(word_count(auxline)))
  call int_alloc(auxline, data)
  if(size(data,1) /= 2) call error('read_msh_dimensions # Wrong number of components')
  pmh%pc(1)%dim = data(2)
  print*, 'Space dimension: '//trim(string(data(2)))

end subroutine


!-----------------------------------------------------------------------
! read_msh_nodes(iu, pmh, line): Read MSH comment section
!-----------------------------------------------------------------------
! iu: msh file unit number
! pmh: pmh mesh
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_nodes(iu, pmh, line)
  integer,                intent(in) :: iu
  character(len=*),       intent(in) :: line
  type(pmh_mesh),      intent(inout) :: pmh ! pmh_mesh
  character(len=MAXPATH)             :: auxline
  character(len=MAXPATH)             :: ch_group,ch_zone,ch_f_indx,ch_l_indx,aux,ch_nd
  integer                            :: zone,f_indx,l_indx
  integer                            :: ios, wc, i, numgroups

  auxline(:) = line(:)
  ! Replace brackets
  call replace_char(auxline, '(', ' ')
  call replace_char(auxline, ')', ' ')
  ! Count words: type may not appear
  wc = word_count(trim(auxline))

  ! Read data as strings
  if(wc == 5) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_nd
    if(ios /= 0) call error('read_msh_nodes # IOerror '//string(ios))
  elseif(wc == 6) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,aux,ch_nd
    if(ios /= 0) call error('read_msh_nodes # IOerror '//string(ios))
  else
    call error('read_msh_nodes # Wrong number of components')
  endif

  ! Converts data from hexadecimal to decimal
  read(ch_zone,fmt='(z20)',iostat=ios) zone
  if(ios /= 0) call error('read_msh_nodes # IOerror '//string(ios))
  if(pmh%pc(1)%dim < int(ch_nd)) pmh%pc(1)%dim = int(ch_nd)
  read(ch_f_indx,fmt='(z20)',iostat=ios) f_indx
  if(ios /= 0) call error('read_msh_nodes # IOerror '//string(ios))
  read(ch_l_indx,fmt='(z20)',iostat=ios) l_indx
  if(ios /= 0) call error('read_msh_nodes # IOerror '//string(ios))

  ! Allocate or extend nodes if needed
  if(.not. allocated(pmh%pc(1)%z)) allocate(pmh%pc(1)%z(pmh%pc(1)%dim,l_indx))
  if(size(pmh%pc(1)%z,1)<pmh%pc(1)%dim .or. size(pmh%pc(1)%z,2)<l_indx) &
    & call extend(pmh%pc(1)%z,pmh%pc(1)%dim,l_indx)

  ! If section contains data
  if(zone /= 0) then
    ! Allocate element group if needed
    if(.not. allocated(pmh%pc(1)%el)) then
      allocate(pmh%pc(1)%el(1))
      numgroups = 1
    else
      call extend_elgroup(pmh%pc(1)%el, size(pmh%pc(1)%el,1)+1)
      numgroups = size(pmh%pc(1)%el,1)
    endif

    if(.not. allocated(pmh%pc(1)%el(numgroups)%ref)) &
      & allocate(pmh%pc(1)%el(numgroups)%ref(l_indx-f_indx+1))
    if(.not. allocated(pmh%pc(1)%el(numgroups)%mm)) &
      & allocate(pmh%pc(1)%el(numgroups)%mm(1,l_indx-f_indx+1))


    pmh%pc(1)%el(numgroups)%type = check_fe(.true., 1, 1, 0, 0)
    pmh%pc(1)%el(numgroups)%nel = 0

    print*, 'Reading zone '//trim(string(zone))//': '//trim(FEDB(pmh%pc(1)%el(numgroups)%type)%desc)

    ! Read node coordinates
    i = f_indx
    do
      if(i > l_indx) exit
      if(is_blank_line(line)) cycle
      read(iu,fmt='(A)', iostat=ios) auxline
      if(ios /= 0) exit
      read(auxline,fmt=*,iostat=ios) pmh%pc(1)%z(:,i)
      pmh%pc(1)%nver = pmh%pc(1)%nver + 1
      pmh%pc(1)%nnod = pmh%pc(1)%nver
      pmh%pc(1)%el(numgroups)%nel = pmh%pc(1)%el(numgroups)%nel + 1
      pmh%pc(1)%el(numgroups)%mm(1,pmh%pc(1)%el(numgroups)%nel) = pmh%pc(1)%nver
      pmh%pc(1)%el(numgroups)%ref(pmh%pc(1)%el(numgroups)%nel) = zone
      i = i + 1
    enddo
  end if

end subroutine


!-----------------------------------------------------------------------
! read_msh_faces(iu, faces, line): Read MSH faces section
!-----------------------------------------------------------------------
! iu: msh file unit number
! faces: array of msh face data
! line: header line of the section
!-----------------------------------------------------------------------
!subroutine read_msh_faces2(iu, faces, line)
!  integer,                         intent(in) :: iu
!  character(len=*),                intent(in) :: line
!  type(msh_faces), allocatable, intent(inout) :: faces(:) ! msh faces
!  character(len=MAXPATH)             :: auxline,ch_group,ch_zone,ch_f_indx
!  character(len=MAXPATH)             :: ch_l_indx,ch_tp,ch_el_tp,ch_aux
!  character(len=MAXPATH), allocatable:: ch_array(:)
!  integer                            :: zone,f_indx,l_indx,tp
!  integer                            :: ios, wc, i, j, numadcells,mix
!
!  auxline(:) = line(:)
!  ! Replace brackets
!  call replace_char(auxline, '(', ' ')
!  call replace_char(auxline, ')', ' ')
!  wc = word_count(trim(auxline))
!
!  ! Count words: type may not appear
!  if(wc == 5) then
!    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_el_tp
!    if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
!  elseif(wc == 6) then
!    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_tp,ch_el_tp
!    if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
!  else
!    call error('read_msh_cells # Wrong number of components')
!  endif
!
!  ! Converts data from hexadecimal to decimal
!  read(ch_zone,fmt='(z20)',iostat=ios) zone
!  if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
!  read(ch_f_indx,fmt='(z20)',iostat=ios) f_indx
!  if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
!  read(ch_l_indx,fmt='(z20)',iostat=ios) l_indx
!  if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
!  tp = int(ch_el_tp)
!
!  ! Allocate or extend faces if needed
!  if(.not. allocated(faces)) then
!    allocate(faces(l_indx))
!  elseif(size(faces) < l_indx) then
!    call extend_faces(faces,l_indx)
!  endif
!
!  ! If section contains data
!  if(zone /= 0) then
!    i = f_indx
!    do
!      if(i > l_indx) exit
!      if(is_blank_line(line)) cycle
!      read(iu,fmt='(A)', iostat=ios) auxline
!      if(ios /= 0) exit
!
!      ! Reads face data into string array
!      wc = word_count(trim(auxline))
!      if(allocated(ch_array)) deallocate(ch_array)
!      allocate(ch_array(wc))
!      read(auxline,fmt=*,iostat=ios) (ch_array(j),j=1,wc)
!
!      ! Read face type if not specified in section header
!      mix = 0
!      if(tp == 0) then; read(ch_array(1),fmt=*,iostat=ios) tp; mix = 1; endif
!      faces(i)%type = get_face_pmh_type(tp)
!      faces(i)%zone = zone
!
!      ! Read face connectivity
!      if(.not. allocated(faces(i)%mm)) allocate(faces(i)%mm(FEDB(faces(i)%type)%lnv))
!      do j=1, FEDB(faces(i)%type)%lnv
!        read(ch_array(j+mix),fmt='(z20)',iostat=ios) faces(i)%mm(j)
!      enddo
!
!      ! Read cells adjacents to a face
!      numadcells = wc - FEDB(faces(i)%type)%lnv - mix
!      if(.not. allocated(faces(i)%adcells)) allocate(faces(i)%adcells(numadcells))
!      do j=1, numadcells
!        read(ch_array(j+mix+FEDB(faces(i)%type)%lnv),fmt='(z20)',iostat=ios) faces(i)%adcells(j)
!      enddo
!
!      i = i + 1
!    enddo
!  endif
!
!end subroutine


!-----------------------------------------------------------------------
! read_msh_faces(iu, faces, line): Read MSH faces section
!-----------------------------------------------------------------------
! iu: msh file unit number
! faces: array of msh face data
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_faces(iu, faces, line)
  integer,                intent(in) :: iu
  character(len=*),       intent(in) :: line
  type(msh_faces),     intent(inout) :: faces ! msh faces
  character(len=MAXPATH)             :: auxline,ch_group,ch_zone,ch_f_indx
  character(len=MAXPATH)             :: ch_l_indx,ch_tp,ch_el_tp
  character(len=MAXPATH), allocatable:: ch_array(:)
  integer                            :: zone,f_indx,l_indx,zone_tp,tp
  integer                            :: ios, wc, i, j, numadcells,mix,ft

  auxline(:) = line(:)
  ! Replace brackets
  call replace_char(auxline, '(', ' ')
  call replace_char(auxline, ')', ' ')
  wc = word_count(trim(auxline))

  ! Count words: type may not appear
  if(wc == 5) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_el_tp
    if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
  elseif(wc == 6) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_tp,ch_el_tp
    if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
  else
    call error('read_msh_faces # Wrong number of components')
  endif

  ! Converts data from hexadecimal to decimal
  read(ch_zone,fmt='(z20)',iostat=ios) zone
  if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
  read(ch_f_indx,fmt='(z20)',iostat=ios) f_indx
  if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
  read(ch_l_indx,fmt='(z20)',iostat=ios) l_indx
  if(ios /= 0) call error('read_msh_faces # IOerror '//string(ios))
  zone_tp = int(ch_el_tp)
  tp = zone_tp

  ! Allocate or extend faces if needed
  if(.not. allocated(faces%type)) then
    allocate(faces%type(l_indx))
    allocate(faces%zone(l_indx))
    allocate(faces%mm(l_indx))
    allocate(faces%adcells(l_indx))
  elseif(size(faces%type) < l_indx) then
    call extend(faces%type,l_indx)
    call extend(faces%zone,l_indx)
    call extend_varlen(faces%mm,l_indx)
    call extend_varlen(faces%adcells,l_indx)
  endif

  ! If section contains data
  if(zone /= 0) then
    if(zone_tp == 0) then
      print*, 'Reading mixed faces zone '//trim(string(zone))
    else
      print*, 'Reading zone '//trim(string(zone))//': '//trim(FEDB(get_face_pmh_type(zone_tp))%desc)
    endif
    i = f_indx
    do
      if(i > l_indx) exit
      if(is_blank_line(line)) cycle
      read(iu,fmt='(A)', iostat=ios) auxline
      if(ios /= 0) exit

      ! Reads face data into string array
      wc = word_count(trim(auxline))
      if(allocated(ch_array)) deallocate(ch_array)
      allocate(ch_array(wc))
      read(auxline,fmt=*,iostat=ios) (ch_array(j),j=1,wc)

      ! Read face type if not specified in section header
      mix = 0
      if(zone_tp == 0) then; read(ch_array(1),fmt=*,iostat=ios) tp; mix = 1; endif
      ft = get_face_pmh_type(tp)

      faces%type(i) = ft
      faces%zone(i) = zone

      ! Read face connectivity
      if(.not. allocated(faces%mm(i)%data)) allocate(faces%mm(i)%data(FEDB(ft)%lnv))
      do j=1, FEDB(ft)%lnv
        read(ch_array(j+mix),fmt='(z20)',iostat=ios) faces%mm(i)%data(j)
      enddo

      ! Read cells adjacents to a face
      numadcells = wc - FEDB(ft)%lnv - mix
      if(.not. allocated(faces%adcells(i)%data)) allocate(faces%adcells(i)%data(numadcells))
      do j=1, numadcells
        read(ch_array(j+mix+FEDB(ft)%lnv),fmt='(z20)',iostat=ios) faces%adcells(i)%data(j)
      enddo

      i = i + 1
    enddo
  endif

end subroutine


!-----------------------------------------------------------------------
! read_msh_cells(iu, cells, line): Read MSH cells section
!-----------------------------------------------------------------------
! iu: msh file unit number
! cells: cell types
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_cells(iu, cells, line)
  integer,                 intent(in) :: iu
  character(len=*),        intent(in) :: line
  type(msh_cells),      intent(inout) :: cells ! cell types
  character(len=MAXPATH)              :: auxline,ch_l_indx,ch_tp,ch_el_tp
  character(len=MAXPATH)              :: ch_group,ch_zone,ch_f_indx
  integer                             :: zone,f_indx,l_indx,tp
  integer                             :: ios, wc, j

  auxline(:) = line(:)
  ! Replace brackets
  call replace_char(auxline, '(', ' ')
  call replace_char(auxline, ')', ' ')

  ! Count words: type may not appear
  wc = word_count(trim(auxline))
  if(wc == 5) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_el_tp
    if(ios /= 0) call error('read_msh_cells # IOerror '//string(ios))
  elseif(wc == 6) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_f_indx,ch_l_indx,ch_tp,ch_el_tp
    if(ios /= 0) call error('read_msh_cells # IOerror '//string(ios))
  else
    call error('read_msh_cells # Wrong number of components')
  endif


  ! Converts data from hexadecimal to decimal
  read(ch_zone,fmt='(z20)',iostat=ios) zone
  if(ios /= 0) call error('read_msh_cells # IOerror '//string(ios))
  read(ch_f_indx,fmt='(z20)',iostat=ios) f_indx
  if(ios /= 0) call error('read_msh_cells # IOerror '//string(ios))
  read(ch_l_indx,fmt='(z20)',iostat=ios) l_indx
  if(ios /= 0) call error('read_msh_cells # IOerror '//string(ios))

  tp = int(ch_el_tp)

  ! Allocate or extend cells if needed
  if(.not. allocated(cells%type)) then
    allocate(cells%type(l_indx))
    allocate(cells%zone(l_indx))
  elseif(size(cells%type) < l_indx) then
    call extend(cells%type,l_indx)
    call extend(cells%zone,l_indx)
  endif

  ! If section contains data
  if(zone /= 0) then
    if(tp == 0) then
      print*, 'Reading mixed cells zone '//trim(string(zone))
    else
      print*, 'Reading zone '//trim(string(zone))//': '//trim(FEDB(get_cell_pmh_type(tp))%desc)
    endif
    ! Homogeneus group of cells
    if(tp /= 0) then
      cells%type(f_indx:l_indx) = tp
    ! Heterogeneus group of cells
    else
        read(iu,fmt=*,iostat=ios) (cells%type(j), j=f_indx, l_indx)
!      i = f_indx
!      cb = 0 ! If exists close brackes exit
!      do
!        if(i > l_indx) exit
!        read(iu,fmt='(A)', iostat=ios) auxline
!        if(ios /= 0) exit
!        if(is_blank_line(auxline)) cycle
!        if(index('0123456', auxline) == 0) exit
!        cb = index(auxline, ')')
!        call replace_char(auxline, '(', ' ')
!        call replace_char(auxline, ')', ' ')
!        wc = word_count(trim(auxline))
!	print*, auxline
!        read(auxline,fmt=*,iostat=ios) (cells%type(j), j=i, i+wx-1)
!        i = i + wc
!        if(cb>0) exit
!      enddo
!      if(i<size(cells%type,1)) call error('reading cell types')
    endif
    cells%zone(f_indx:l_indx) = zone
  endif

end subroutine

!-----------------------------------------------------------------------
! read_msh_zones(iu, zones, line): Read MSH zone section
!-----------------------------------------------------------------------
! iu: msh file unit number
! zones: array of zone info
! line: header line of the section
!-----------------------------------------------------------------------
subroutine read_msh_zones(izones, line)
  character(len=*),        intent(in) :: line
  type(msh_zone),       intent(inout) :: izones ! interior zones
  character(len=MAXPATH)              :: ch_group,ch_zone,ch_type,ch_name
  character(len=MAXPATH)              :: auxline
  integer                             :: ios, wc, i

  auxline(:) = line(:)
  ! Replace brackets
  call replace_char(auxline, '(', ' ')
  call replace_char(auxline, ')', ' ')

  ! Count words: type may not appear
  wc = word_count(trim(auxline))
  if(wc == 4) then
    read(auxline,fmt=*, iostat=ios) ch_group,ch_zone,ch_type,ch_name
    if(ios /= 0) call error('read_msh_cells # IOerror '//string(ios))
  else
    call error('read_msh_cells # Wrong number of components')
  endif

!  if(trim(ch_type) == 'interior') then
    if(.not. allocated(izones%id)) then
      i = 1
      allocate(izones%id(i))
      allocate(izones%names(i))
      allocate(izones%types(i))
    else
      i = size(izones%id,1)+1
      call extend(izones%id,i)
      call extend(izones%names,i)
      call extend(izones%types,i)
    endif

    izones%id(i) =  int(ch_zone)
    izones%names(i) =  trim(ch_name)
    izones%types(i) =  trim(ch_type)
!  endif

end subroutine


end module
