module module_xread
!-----------------------------------------------------------------------
! Module to read marked data files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 04/03/2012
!
! PUBLIC PROCEDURES:
!   xread: reads scalar data between marks
!   xread_alloc: reads array data between marks
!   xlist: lists all the marks enclosed in a given path
!   is_alphanumeric: checks whether the first character is alphanumeric
!
! EXAMPLE:
!   integer :: i
!   logical :: res
!   res = xread(10, 'leaf', i)
!   !A valid xml tag would be: <leaf> 23 </leaf>
!
! NOTE:
!   Numerical values must be separated between them and from the marks that 
!   enclose them by spaces or returns; character values must be separated 
!   between them and from the marks that enclose them only by returns.
!
!   When marks starts by an alphanumeric character, xread follows the
!   path of tags until the last mark; otherwise the mark is searched 
!   independently of the xml hierarchy  
!-----------------------------------------------------------------------
use module_compiler_dependant, only: iostat_end, real64
use module_report, only: info, error
use module_convers, only: maxpath, string, adjustlt, lcase, int_count, dble_count, word
implicit none

!Private variables
integer, private :: nlines !counter for lines already read

!Private procedures
private :: xread_int, xread_real64, xread_char
private :: xread_int_alloc, xread_real64_alloc, xread_char_alloc
private :: search_once, search_mark, search_path, search_close_mark
private :: last_mark

!Interfaces
interface xread; module procedure xread_int;          end interface
interface xread; module procedure xread_real64;       end interface
interface xread; module procedure xread_char;         end interface
interface xread; module procedure xread_int_alloc;    end interface
interface xread; module procedure xread_real64_alloc; end interface
interface xread; module procedure xread_char_alloc;   end interface
 
contains

!-----------------------------------------------------------------------
! xread: read an integer scalar
!-----------------------------------------------------------------------
function xread_int(id, mark, d, rew) result(res)
integer, intent(in) :: id
character(*), intent(in) :: mark
integer, intent(out) :: d
logical, optional, intent(in) :: rew
character(maxpath) :: str
integer :: res

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
backspace(id)
read (id, *, iostat=res) str, d
if (res /= 0) call error('(xread_int/read) error while reading mark "'//&
trim(mark)//'", #'//trim(string(res)))
call info('read <'//trim(adjustlt(mark))//'>: '//string(d))
end function

!-----------------------------------------------------------------------
! xread: read a real64 scalar
!-----------------------------------------------------------------------
function xread_real64(id, mark, d, rew) result(res)
integer, intent(in) :: id
character(*), intent(in) :: mark
real(real64), intent(out) :: d
logical, optional, intent(in) :: rew
character(maxpath) :: str
integer :: res

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
backspace(id)
read (id, *, iostat=res) str, d
if (res /= 0) call error('(xread_real64/read) error while reading mark "'//&
trim(mark)//'", #'//trim(string(res))) 
call info('read <'//trim(adjustlt(mark))//'>: '//string(d))
end function

!-----------------------------------------------------------------------
! xread: read a character scalar
!-----------------------------------------------------------------------
function xread_char(id, mark, d, rew) result(res)
integer, intent(in) :: id
character(*), intent(in) :: mark
character(*), intent(out) :: d
logical, optional, intent(in) :: rew
integer :: res

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
read (id, '(a)', iostat=res) d
if (res /= 0) call error('(xread_char/read) error while reading mark "'//&
trim(mark)//'", #'//trim(string(res))) 
call info('read <'//trim(adjustlt(mark))//'>: '//'"'//trim(d)//'"')
end function

!-----------------------------------------------------------------------
! xread: read an integer array
!-----------------------------------------------------------------------
function xread_int_alloc(id, mark, d, rew) result(res)
integer, intent(in) :: id
character(*), intent(in) :: mark
integer, intent(out), allocatable :: d(:)
logical, optional, intent(in) :: rew
character(maxpath) :: str, cad
integer :: res, n, i

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
backspace(id)
!read records until ending mark and count words
n = 0; nlines = 0
do
  read (id, '(a)', iostat=res) str
  nlines = nlines + 1
  if (res /= 0) call error('(xread_int_alloc/read) error while reading mark "'//&
  trim(mark)//'", #'//trim(string(res))) !error found
  n = n + int_count(str)
  if (index(lcase(str), lcase('</'//trim(adjustlt(last_mark(mark)))//' ')) > 0 .or. &
      index(lcase(str), lcase('</'//trim(adjustlt(last_mark(mark)))//'>')) > 0) exit !end mark found
end do
do i = 1, nlines; backspace(id); end do
!allocation
if (n <= 0) call error('(xread_int_alloc/read) no data was read in mark "'//trim(mark)//'".')
allocate(d(n), stat = res, errmsg = cad)
if (res /= 0) call error('(xread_int_alloc) unable to allocate variable: '//trim(cad))
!read records
read (id, *, iostat=res) str, d
if (res /= 0) call error('(xread_int_alloc/read) error while reading content of mark "'//&
trim(mark)//'", #'//trim(string(res))) !error found
call info('read <'//trim(adjustlt(mark))//'>: '//trim(string(d)))
end function

!-----------------------------------------------------------------------
! xread: read a real64 array
!-----------------------------------------------------------------------
function xread_real64_alloc(id, mark, d, rew) result(res)
integer, intent(in) :: id
character(*), intent(in) :: mark
real(real64), intent(out), allocatable :: d(:)
logical, optional, intent(in) :: rew
character(maxpath) :: str, cad
integer :: res, n , i

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
backspace(id)
!read records until ending mark and count words
n = 0; nlines = 0
do
  read (id, '(a)', iostat=res) str
  nlines = nlines + 1
  if (res /= 0) call error('(xread_real64_alloc/read) error while reading mark "'//&
  trim(mark)//'", #'//trim(string(res))) !error found
  n = n + dble_count(str)
  if (index(lcase(str), lcase('</'//trim(adjustlt(last_mark(mark)))//' ')) > 0 .or. &
      index(lcase(str), lcase('</'//trim(adjustlt(last_mark(mark)))//'>')) > 0) exit !end mark found
end do
do i = 1, nlines; backspace(id); end do
!allocation
if (n <= 0) call error('(xread_real64_alloc/read) no data was read in mark "'//trim(mark)//'".')
allocate(d(n), stat = res, errmsg = cad)
if (res /= 0) call error('(xread_real64_alloc) unable to allocate variable: '//trim(cad))
!read records
read (id, *, iostat=res) str, d
if (res /= 0) call error('(xread_real64_alloc/read) error while reading content of mark "'//&
trim(mark)//'", #'//trim(string(res))) !error found
call info('read <'//trim(adjustlt(mark))//'>: '//trim(string(d)))
end function

!-----------------------------------------------------------------------
! xread: read a character array
!-----------------------------------------------------------------------
function xread_char_alloc(id, mark, d, rew) result(res)
integer, intent(in) :: id
character(*), intent(in) :: mark
character(*), intent(out), allocatable :: d(:)
logical, optional, intent(in) :: rew
character(maxpath) :: str, cad
character(4096) :: longstr
integer :: res, n, i

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
!read records until ending mark and count records
n = 0; nlines = 0
do
  read (id, '(a)', iostat=res) str
  nlines = nlines + 1
  if (res /= 0) call error('(xread_char_alloc/read) error while reading mark "'//&
  trim(mark)//'", #'//trim(string(res))) !error found
  if (index(lcase(str), lcase('</'//trim(adjustlt(last_mark(mark)))//' ')) > 0 .or. &
      index(lcase(str), lcase('</'//trim(adjustlt(last_mark(mark)))//'>')) > 0) exit !end mark found
  n = n + 1
end do
do i = 1, nlines; backspace(id); end do
!allocation
if (n <= 0) call error('(xread_char_alloc/read) no data was read in mark "'//trim(mark)//'".')
allocate(d(n), stat = res, errmsg = cad)
if (res /= 0) call error('xread_char_alloc) unable to allocate variable: '//trim(cad))
!read records
do i = 1, n
  read (id, '(a)', iostat=res) d(i)
  if (res /= 0) call error('(xread_char_alloc/read) error while reading '//trim(string(i))//'th element of mark "'//&
  trim(mark)//'", #'//trim(string(res))) !error found
end do
!print result
longstr = '"'//trim(d(1))//'"'
do i = 2, n
  longstr = trim(longstr)//', "'//trim(d(i))//'"'
end do  
call info('read <'//trim(adjustlt(mark))//'>: '//trim(longstr))
end function

!-----------------------------------------------------------------------
! xlist: lists all the marks enclosed in a given path
!-----------------------------------------------------------------------
subroutine xlist(id, path, val, rew)
integer,                   intent(in)  :: id
character(*),              intent(in)  :: path
character(*), allocatable, intent(out) :: val(:)
logical, optional,         intent(in)  :: rew
character(maxpath) :: mark, str
character(4096) :: longstr
integer :: n, res, i

!find the last mark
mark = last_mark(path)
!follow the path
res = search_path(id, path, rew)
n = 0; nlines = 0
do
  !search marks inside
  read (id, '(a)', iostat=res) str; str = adjustl(str)
  nlines = nlines + 1
  if (res == iostat_end) exit !end-of-file found
  if (res /= 0) call error('(xlist/read) error while reading mark "'//trim(mark)//'", #'//&
  &trim(string(res))) !error found
  if (index(lcase(str), '</'//trim(mark)//' ') == 1 .or. &
      index(lcase(str), '</'//trim(mark)//'>') == 1) exit !close mark found
  if (index(str, '<') == 1) then !member-mark found
    backspace(id); nlines = nlines - 1 
    res = search_close_mark(id, aword(str(2:len_trim(str))))
    if (res == iostat_end) return !end-of-file before the close member-mark was found
    n = n + 1 
  end if
end do
do i = 1, nlines; backspace(id); end do
if (res /= 0 .or. n <= 0) return
!write member-marks
allocate(val(n), stat = res, errmsg = str)
if (res /= 0) call error('(xlist) unable to allocate variable: '//trim(str))
n = 0; nlines = 0
do while (n < size(val,1))
  !search marks inside
  read (id, '(a)', iostat=res) str; str = adjustlt(str)
  nlines = nlines + 1
  if (index(str, '<') == 1) then !member-mark found
    n = n + 1 
    val(n) = aword(str(2:len_trim(str)))
    backspace(id); nlines = nlines - 1 
    res = search_close_mark(id, val(n))
  end if
end do
do i = 1, nlines; backspace(id); end do
!print result
if (n <= 0) then
  call info('list <'//trim(adjustlt(path))//'>: no members inside')
else
  longstr = '"'//trim(val(1))//'"'
  do i = 2, n
    longstr = trim(longstr)//', "'//trim(val(i))//'"'
  end do  
  call info('list <'//trim(adjustlt(path))//'>: '//trim(longstr))
end if
end subroutine

!-----------------------------------------------------------------------
! is_alphanumeric: checks whether the first character is alphanumeric
!-----------------------------------------------------------------------
function is_alphanumeric(c) result(res)
character(*) :: c
logical :: res

if ((iachar('a') <= iachar(c) .and. iachar(c) <= iachar('z')) .or. &
    (iachar('A') <= iachar(c) .and. iachar(c) <= iachar('Z')) .or. &
    (iachar('0') <= iachar(c) .and. iachar(c) <= iachar('9'))) then
  res = .true.
else
  res = .false.
end if
end function

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! search_once: searches for a mark only once
! RETURN: 0 if the mark is found; iostat_end if end-of-file is found; 
! non-zero otherwise
!-----------------------------------------------------------------------
function search_once(id, mark) result(res)
integer,      intent(in) :: id
character(*), intent(in) :: mark
character(maxpath) :: str
integer :: res

read (id, '(a)', iostat=res) str
nlines = nlines + 1
if (res == iostat_end) return !end-of-file found
if (res /= 0) call error('(search_mark_once/read) error while reading mark "'//trim(mark)//'", #'//&
trim(string(res))) !error found
if (index(adjustlt(lcase(str)), trim(adjustlt(lcase(mark)))//' ') /= 1 .and. &
    index(adjustlt(lcase(str)), trim(adjustlt(lcase(mark)))//'>') /= 1 ) res = 1 !mark not found
end function

!-----------------------------------------------------------------------
! search_mark: searches for '<mark' in the unit id
! RETURN: 0 if the mark is found; iostat_end if end-of-file is found; 
! non-zero otherwise
!-----------------------------------------------------------------------
function search_mark(id, mark) result(res)
integer,      intent(in) :: id
character(*), intent(in) :: mark
integer :: res

do
  res = search_once(id, '<'//trim(adjustlt(lcase(mark))))
  if (res == iostat_end .or. res == 0) return !end-of-file or mark found
end do
end function

!-----------------------------------------------------------------------
! search_path: searchs a mark or a path; if mark(1:1) is alphanumeric, 
! a mark is searched with procedure search_mark; otherwise, 'mark' is
! considered a path of marks separated by mark(1:1)
! RETURN: 0 if the mark is found; iostat_end if end-of-file is found; 
! non-zero otherwise
!-----------------------------------------------------------------------
function search_path(id, mark, rew) result(res)
integer,      intent(in) :: id
character(*), intent(in) :: mark
logical, optional, intent(in) :: rew
character(maxpath) :: part, path
character(1) :: sep
integer :: res, p

!optional rewind
if (present(rew)) then
  if (rew) rewind(id)
end if
if (is_alphanumeric(mark)) then
  !a single mark is provided
  res = search_mark(id, mark)
else
  !a path is provided
  sep = mark(1:1) !separator
  path = mark(2:len_trim(mark)) !it will contain the remaining path
  part = path     !it will contain the first mark of the remaining part
  do while (len_trim(part) > 0)
    p = index(path, sep)
    if (p > 0) then
      part = path(:p-1); path = path(p+1:len_trim(path))
    else
      part = path;       path = ' '
    endif
    res = search_mark(id, part)
    part = path
  enddo
end if

end function

!-----------------------------------------------------------------------
! search_close_mark: search a close mark
!-----------------------------------------------------------------------
function search_close_mark(id, mark) result(res)
integer,      intent(in) :: id
character(*), intent(in) :: mark
character(maxpath) :: str
integer :: res, n

n = 0 !number of open marks
do 
  read (id, '(a)', iostat=res) str
  nlines = nlines + 1
  if (res == iostat_end) return !end-of-file found
  if (res /= 0) call error('(search_close_mark/read) error while reading mark "'//trim(mark)//'", #'//&
  &trim(string(res))) !error found
  if (index(adjustlt(lcase(str)),  '<'//trim(adjustlt(lcase(mark)))//' ') == 1 .or. &
      index(adjustlt(lcase(str)),  '<'//trim(adjustlt(lcase(mark)))//'>') == 1 ) n = n + 1 !open  mark found
  if (index(adjustlt(lcase(str)), '</'//trim(adjustlt(lcase(mark)))//' ') > 0 .or. &
      index(adjustlt(lcase(str)), '</'//trim(adjustlt(lcase(mark)))//'>') > 0 )  n = n - 1 !close mark found
  if (n == 0) exit
enddo
end function

!-----------------------------------------------------------------------
! last_mark: find the last mark of a path
!-----------------------------------------------------------------------
function last_mark(path) result(mark)
character(*),  intent(in) :: path
character(len(path)) :: mark
integer :: p

if (is_alphanumeric(path)) then
  mark = trim(adjustlt(lcase(path)))
else
  p = index(path, path(1:1), back=.true.)
  mark = trim(adjustlt(lcase(path(p+1:len_trim(path)))))
end if
end function

!--------------------------------------------------------------------
! word: gets the first word of a string removing non alphanumeric 
! characters at the end
! NOTE: please use trim(word(cad))
!--------------------------------------------------------------------
function aword(cad) result(res)
character(*), intent(in) :: cad
character(maxpath) :: res
integer :: ios
  
read(cad, *, iostat=ios) res
if (ios /= 0) call error('aword/read, #'//trim(string(ios)))
res = adjustl(res)
do while (.not. is_alphanumeric(res(len_trim(res):len_trim(res))))
  res = res(1:len_trim(res)-1)
end do  
end function

end module
