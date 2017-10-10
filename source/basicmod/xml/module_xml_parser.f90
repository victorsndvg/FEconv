module module_xml_parser
!-----------------------------------------------------------------------
! Module for reading xml files
! Last update: 26/04/2009
!-----------------------------------------------------------------------
use module_ALLOC
use module_CONVERS
use module_REPORT
use module_FILES
implicit none

!Constants
character(len=*), dimension(4), parameter, private ::  OPEN_MARK_MEMBERS = &
(/ '<menu   ', '<submenu','<struct ', '<leaf   ' /)
character(len=*), dimension(4), parameter, private :: CLOSE_MARK_MEMBERS = &
(/ '</menu   ', '</submenu','</struct ', '</leaf   ' /)
character(len=*), dimension(2), parameter, private ::  OPEN_MARK_LEAF = (/ '<leaf', '<data'/)
character(len=*), dimension(2), parameter, private :: CLOSE_MARK_LEAF = (/ '</leaf', '</data'/)
character(len=*), dimension(1), parameter, private ::  OPEN_MARK_ELEMENT = (/ '<element' /)
character(len=*), dimension(1), parameter, private :: CLOSE_MARK_ELEMENT = (/ '</element' /)

!Class attributes

!Private methods
private :: search_mark_once, search_mark, search_close_mark, follow_path, &
           last_part, cut_end_delimiter
private :: fread_real, fread_vreal, fread_vreal_alloc, &
           fread_complex, fread_vcomplex, fread_vcomplex_alloc, &
           fread_char, fread_vchar, fread_vchar_alloc

!Interfaces
interface fread; module procedure fread_real;  end interface
interface fread; module procedure fread_vreal; end interface
interface fread; module procedure fread_complex;  end interface
interface fread; module procedure fread_vcomplex; end interface
interface fread; module procedure fread_char;  end interface
interface fread; module procedure fread_vchar; end interface
interface fread_alloc; module procedure fread_vreal_alloc; end interface
interface fread_alloc; module procedure fread_vcomplex_alloc; end interface
interface fread_alloc; module procedure fread_vchar_alloc; end interface

contains

!-----------------------------------------------------------------------
! fopen: open a xml file
!-----------------------------------------------------------------------
function fopen(datxml) result(un)

character(len=*), intent(in), optional :: datxml !file name
character(len=MAXPATH) :: xmlfile, arg1
integer :: un !associated unit number
integer :: ios, length, status

!find a valid xmlfile
xmlfile = ' '
if (present(datxml)) then
  xmlfile = datxml
else
  if (command_argument_count() == 2) then
    call get_command_argument(1, arg1, length, status)
    if (status /= 0) call error('fopen/get_command_argument, '&
     &//'the first command argument cannot be read')
    if (trim(adjustl(arg1)) == '-xml') then
      call get_command_argument(2, xmlfile, length, status)
      if (status /= 0) call error('fopen/get_command_argument, '&
       &//'the second command argument cannot be read')
    else
      call error('fopen/get_command_argument, '&
       &//'the first command argument is not recognized (must be -xml)')
    endif
  elseif (command_argument_count() == 1) then
    call get_command_argument(1, arg1, length, status)
    if (status /= 0) call error('fopen/get_command_argument, '&
     &//'the first command argument cannot be read')
    if (trim(adjustl(arg1)) == '-xml') then
      xmlfile = 'local.dat.xml'
    else
      call error('fopen/get_command_argument, '&
       &//'the first command argument is not recognized (must be -xml)')
    endif
  elseif (command_argument_count() == 0) then
    xmlfile = 'local.dat.xml'
  else
    call error('fopen/get_command_argument, '&
     &//'too many arguments')
  end if
end if
if (len_trim(xmlfile)==0) call error('fopen, filename is empty')
call info('fopen (xmlfile), '//trim(xmlfile))
! open xmlfile
un = get_unit()
open (unit=un, file=xmlfile, form='formatted', iostat=ios, &
status='old', position='rewind')
if (ios /= 0) call error('fopen, #'//trim(string(ios)))

end function

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_real(un, path, var)

integer,          intent(in)  :: un
character(len=*), intent(in)  :: path
real(real64),     intent(out) :: var
integer :: res, tn = 1

call follow_path(un, path, back = .true.)
if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='float', tnum=tn)) call error(trim(path)//'), not found')
!get_elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
read (un, *, iostat=res) var
if (res /= 0) call error('fread_real/read ('//trim(path)//'), #'//trim(string(res)))
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
call info('fread_real ('//trim(path)//'), '//trim(string(var)))

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vreal(un, path, var)

integer,          intent(in)  :: un
character(len=*),           intent(in)  :: path
real(real64), dimension(:), intent(out) :: var
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='float', tnum=tn)) call error(trim(path)//'), not found')
if (tn > size(var,1)) call error('fread_vreal ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vreal/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  call info('fread_vreal ('//trim(path)//'), '//trim(string(var(i))))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vreal_alloc(un, path, var, realloc)

integer,                    intent(in)  :: un
character(len=*),           intent(in)  :: path
real(real64), dimension(:), intent(inout), allocatable :: var
logical,                    intent(in), optional :: realloc
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='float', tnum=tn)) call error(trim(path)//'), not found')
if (present(realloc)) then; if (realloc) then
  if (allocated(var)) call dealloc(var)
endif; endif
if (.not. allocated(var)) call alloc(var, tn)
if (tn > size(var,1)) call error('fread_vreal_alloc ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vreal/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  call info('fread_vreal_alloc ('//trim(path)//'), '//trim(string(var(i))))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_complex(un, path, var)

integer,          intent(in)  :: un
character(len=*), intent(in)  :: path
character(len=128)            :: tempstring
complex(real64),  intent(out) :: var
integer :: res, tn = 1

call follow_path(un, path, back = .true.)
if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='complex', tnum=tn)) call error(trim(path)//'), not found')
!get_elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
read (un, *, iostat=res) var
if (res /= 0) call error('fread_complex/read ('//trim(path)//'), #'//trim(string(res)))
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
write(tempstring,*) var
call info('fread_complex ('//trim(path)//'), '//trim(tempstring))

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vcomplex(un, path, var)

integer,          intent(in)  :: un
character(len=*),           intent(in)  :: path
character(len=128)                      :: tempstring
complex(real64), dimension(:), intent(out) :: var
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='complex', tnum=tn)) call error(trim(path)//'), not found')
if (tn > size(var,1)) call error('fread_vcomplex ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vcomplex/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  write(tempstring,*) var(i)
  call info('fread_vcomplex ('//trim(path)//'), '//trim(tempstring))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vcomplex_alloc(un, path, var, realloc)

integer,                    intent(in)  :: un
character(len=*),           intent(in)  :: path
character(len=128)                      :: tempstring
complex(real64), dimension(:), intent(inout), allocatable :: var
logical,                    intent(in), optional :: realloc
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='complex', tnum=tn)) call error(trim(path)//'), not found')
if (present(realloc)) then; if (realloc) then
  if (allocated(var)) deallocate(var)
endif; endif
if (.not. allocated(var)) allocate(var(tn))
if (tn > size(var,1)) call error('fread_vcomplex_alloc ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vcomplex_alloc/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  write(tempstring,*) var(i)
  call info('fread_vcomplex_alloc ('//trim(path)//'), '//trim(tempstring))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_char(un, path, var)

integer,          intent(in)  :: un
character(len=*), intent(in)  :: path
character(len=*), intent(out) :: var
integer :: res, tn = 1

call follow_path(un, path, back = .true.)
if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), tnum=tn)) & !typ can be diverse
  call error(trim(path)//'), not found')
!get_elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
read (un, '(a)', iostat=res) var
var = adjustlt(var)
if (res /= 0) call error('fread_char/read ('//trim(path)//'), #'//trim(string(res)))
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
call info('fread_char ('//trim(path)//'), '//trim(var))

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vchar(un, path, var)

integer,                        intent(in)    :: un
character(len=*),               intent(in)    :: path
character(len=*), dimension(:), intent(inout) :: var
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not. search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), tnum=tn)) & !typ can be diverse
  call error('fread_vchar ('//trim(path)//'), not found')
if (tn > size(var,1)) call error('fread_vchar ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, '(a'//trim(string(len(var(1))))//')', iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vreal/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  var(i) = adjustlt(var(i))
  call info('fread_vchar ('//trim(path)//'), '//trim(var(i)))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vchar_alloc(un, path, var, realloc)

integer,          intent(in)  :: un
character(len=*), intent(in) :: path
character(len=*), intent(inout), dimension(:), allocatable :: var
logical,          intent(in), optional :: realloc
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not. search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), tnum=tn)) & !typ can be diverse
  call error('fread_vchar_alloc ('//trim(path)//'), not found')
if (present(realloc)) then; if (realloc) then
  if (allocated(var)) call dealloc(var)
endif; endif
if (.not. allocated(var)) call alloc(var, tn)
if (tn > size(var,1)) call error('fread_vreal_alloc ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, '(a'//trim(string(len(var(1))))//')', iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vchar/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  var(i) = adjustlt(var(i))
  call info('fread_vchar_alloc ('//trim(path)//'), '//trim(var(i)))
enddo

end subroutine

!-----------------------------------------------------------------------
! flist: list the members of a structure
!-----------------------------------------------------------------------
subroutine flist(un, path, var)

integer,               intent(in)  :: un
character(len=*),      intent(in)  :: path
character(len=*), allocatable, dimension(:) :: var
character(len=10*MAXPATH), allocatable, dimension(:) :: tmp
integer :: p, n, res
character(len=10*MAXPATH) :: line
logical :: close_mark_found

n = 0
!initial allocation of tmp
if (.not. allocated(tmp)) call alloc(tmp, 10)
!follow the path (advancing the last line)
call follow_path(un, path)
!check whether a close mark is reached
close_mark_found = search_mark_once(un, path, CLOSE_MARK_MEMBERS, back = .true.)
!loop to catch members
do while (.not. close_mark_found)
  !searchs a new member
  call search_mark(un, path, OPEN_MARK_MEMBERS, back=.true.)
  !stores the name
  read (un, fmt='(A)', iostat=res) line
  if (res /= 0) call error('flist/read ('//trim(path)//'), #'//trim(string(res)))
  p = index(line, 'name=')
  n = n + 1; if (n > size(tmp, 1)) call extend(tmp, n, fit=.false.)
  tmp(n) = trim(string(cut_end_delimiter(line(p+5:),'>')))
  !searchs the end of the member
  call search_close_mark(un, path)
  !check whether a close mark is reached
  close_mark_found = search_mark_once(un, path, CLOSE_MARK_MEMBERS, back = .true.)
enddo
!advance one line in file
close_mark_found = search_mark_once(un, path, CLOSE_MARK_MEMBERS)
!final result
if (allocated(var)) call dealloc(var)
call alloc(var, n)
var(1:n) = tmp(1:n)

end subroutine

!-----------------------------------------------------------------------
! fclose: close a xml file
!-----------------------------------------------------------------------
subroutine fclose(un)

integer, intent(in)  :: un
integer :: ios
close(unit=un, iostat=ios)
if (ios /= 0) call error('read_xml/close, #'//trim(string(ios)))

end subroutine

!***********************************************************************
! PRIVATE  PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! search_mark_once: searchs a mark only once
! RETURN: .true. if the mark is found
!         .false. otherwise
!-----------------------------------------------------------------------
recursive function search_mark_once(un, path, marks, name, typ, tnum, advance, back) result(res)
integer,          intent(in) :: un
character(len=*), intent(in) :: path
character(len=*), intent(in), dimension(:) :: marks
character(len=*), intent(in), optional :: name, typ
integer,          intent(inout), optional :: tnum
logical, intent(in), optional :: back, advance
logical :: res
integer :: ios, i, p
character(len=10*MAXPATH) :: line

res = .false.
!read a line
read (un, fmt='(A)', iostat=ios) line
if (ios /= 0) call error('search_mark_once/read ('//trim(path)//'), #'//trim(string(ios)))
!backspace
if (present(back)) then
  if (back) then
    backspace(unit=un, iostat=ios)
    if (ios /= 0) call error('search_mark_once/backspace ('//trim(path)//'), #'//trim(string(ios)))
  endif
endif
do i = 1, size(marks, 1)
  if (index(line, trim(marks(i))) > 0) then
    !check name
    if (present(name)) then
      p = index(line, 'name=')
      if (trim(name) /= trim(string(cut_end_delimiter(line(p+5:),'>')))) then
        if (present(advance)) then
          if (advance) call search_close_mark(un, path)
        endif  
        cycle
      endif
    endif
    !check type
    if (present(typ)) then
      p = index(line, 'type=')
      if (trim(typ) /= trim(string(cut_end_delimiter(line(p+5:),'>')))) then
        if (present(advance)) then
          if (advance) call search_close_mark(un, path)
        endif  
        cycle
      endif
    endif
    !check totalnum
    if (present(tnum)) then
      p = index(line, 'totalnum=')
      if (tnum > 0) then
        if (tnum /= int(string(cut_end_delimiter(line(p+9:),'>')))) then
          if (present(advance)) then
            if (advance) call search_close_mark(un, path)
          endif  
          cycle
        endif
      else
        tnum = int(string(cut_end_delimiter(line(p+9:),'>')))
      endif
    endif
    !the mark, name, type and/or tnum matches
    res = .true.; return
  endif
enddo

end function

!-----------------------------------------------------------------------
! search_mark: searchs a mark
!-----------------------------------------------------------------------
subroutine search_mark(un, path, marks, name, typ, tnum, back, advance)
integer,          intent(in)  :: un
character(len=*), intent(in) :: path
character(len=*), intent(in), dimension(:) :: marks
character(len=*), intent(in), optional :: name, typ
integer,          intent(inout), optional :: tnum
logical, intent(in), optional :: back, advance
integer :: ios

do
  if (search_mark_once(un, path, marks, name, typ, tnum, advance, back=.false.)) then
    !backspace
    if (present(back)) then; if (back) then
      backspace(unit=un, iostat=ios)
      if (ios /= 0) call error('search_mark/backspace, ('//trim(path)//') #'//trim(string(ios)))
    endif; endif
    return
  endif
  !ends the loop if a close mark is found
  if (search_mark_once(un, path, CLOSE_MARK_MEMBERS, back = .true.)) exit
enddo
!mark not found
call error('search_mark ('//trim(path)//'), not found')

end subroutine

!-----------------------------------------------------------------------
! search_close_mark: search a close mark
!-----------------------------------------------------------------------
recursive subroutine search_close_mark(un, path)

integer,           intent(in)  :: un
character(len=*),  intent(in) :: path
integer :: n

n = 1 !number of open marks
do while (n > 0)
  if (search_mark_once(un, path, OPEN_MARK_MEMBERS, back=.true.)) n = n + 1
  if (search_mark_once(un, path, CLOSE_MARK_MEMBERS)) n = n - 1
enddo

end subroutine

!-----------------------------------------------------------------------
! follow_path: follow the path
!-----------------------------------------------------------------------
subroutine follow_path(un, path, back)

integer,           intent(in)  :: un
character(len=*),  intent(in) :: path
logical, optional, intent(in) :: back
character(len=len(path)) :: lpath, parte
character(len=1) :: separador
integer :: p, res

rewind(un)
separador = path(1:1)
lpath = path(2:)
parte = lpath
do while (len_trim(parte) > 0)
  p = index(lpath, separador)
  if (p > 0) then
    parte = lpath(:p-1); lpath = lpath(p+1:)
  else
    parte = lpath;       lpath = ' '
  endif
  call search_mark(un, path, OPEN_MARK_MEMBERS, parte, advance=.true.)
  parte = lpath
enddo
if (present(back)) then; if (back) then
  backspace(unit=un, iostat=res)
  if (res /= 0) call error('follow_path/backspace ('//trim(path)//'), #'//trim(string(res)))
endif; endif

end subroutine

!-----------------------------------------------------------------------
! last_part: extracts the last part of a path
!-----------------------------------------------------------------------
function last_part(path) result(res)

character(len=*),  intent(in) :: path
character(len=len(path)) :: res
character(len=1) :: separador
integer :: p

separador = path(1:1)
p = index(path, separador, back=.true.)
res = path(p+1:)

end function

!-----------------------------------------------------------------------
! cut_end_delimiter: cuts an end delimiter
!-----------------------------------------------------------------------
function cut_end_delimiter(str, delimiter) result(res)
character(len=*), intent(in) :: str, delimiter
character(len=len(str)) :: res
integer :: p

res = str
p = index(str, delimiter)
if  (p > 0) res = str(:p-1)

end function

end module

