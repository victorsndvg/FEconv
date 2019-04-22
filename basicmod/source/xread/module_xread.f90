module module_xread_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to read simple markup language.  
!!
!! Files can have any kind of marks composed of alphanumeric characters. Elements of numeric arrays must be separated, 
!! among them and from the marks, either with blank spaces or returns. Elements of character arrays must be separated, 
!! among them and from the marks, with returns.  
!!
!! __Example:__ _data.xml_  
!! `<?xml version="1.0" encoding="UTF-8" ?>`  
!! `<data>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<cosa> 2 </cosa>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<dato> 1 2 3 4 5 6 </dato>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<constant> 0.5 </constant>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<c> 0.5 0.7 </c>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<mesh>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`a mesh`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</mesh>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<m>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`a mesh`  
!! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`somethig else`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</m>`  
!! `</data>`  
!
! PUBLIC PROCEDURES:
!   xread: reads scalar data between marks
!   xread_alloc: reads array data between marks
!   xlist: lists all the marks enclosed in a given path
!   is_alphanumeric: checks whether the first character is alphanumeric
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: iostat_end, real64
use module_report_bmod, only: info, error
use module_convers_bmod, only: maxpath, string, adjustlt, lcase, int_count, dble_count, word
implicit none

!Private variables
integer, private :: nlines !! Counter for already-read lines.

!Private procedures
private :: xread_int, xread_real64, xread_char
private :: xread_int_alloc, xread_real64_alloc, xread_char_alloc
private :: search_once, search_mark, search_close_mark
private :: last_mark

!Interfaces
interface xread
  !! Read simple markup language.  
  !! __Example:__  
  !! `program test`  
  !! `use basicmod, only: maxpath, xread, xlist`  
  !! `implicit none`  
  !! `integer :: res, a`  
  !! `integer, allocatable :: i(:)`  
  !! `real(real64) :: x`  
  !! `real(real64), allocatable :: x2(:)`  
  !! `character(maxpath) :: c`  
  !! `character(maxpath), allocatable :: c2(:), val(:)`  
  !! `open(10, file='data.xml')`  
  !! `!  `  
  !! `! search unique marks`  
  !! `res = xread(10, 'cosa',     a,  rew=.true.)`  
  !! `res = xread(10, 'dato',     i,  rew=.true.)`  
  !! `res = xread(10, 'constant', x,  rew=.true.)`  
  !! `res = xread(10, 'c',        x2, rew=.true.)`  
  !! `res = xread(10, 'mesh',     c,  rew=.true.)`  
  !! `res = xread(10, 'm',        c2, rew=.true.)`  
  !! `close(10)`  
  !! `!  `  
  !! `! search mark paths`  
  !! `open(10, file='data.xml')`  
  !! `res = xread(10, '/data/cosa',     a,  rew=.true.)`  
  !! `res = xread(10, '/data/dato,      i,  rew=.true.)`  
  !! `res = xread(10, '/data/constant', x,  rew=.true.)`  
  !! `res = xread(10, '/data/c',        x2, rew=.true.)`  
  !! `res = xread(10, '/data/mesh',     c,  rew=.true.)`  
  !! `res = xread(10, '/data/m',        c2, rew=.true.)`  
  !! `close(10)`  
  !! `end program`  
  !! @note If the first character of `mark` is alphanumeric, it is a separator and `mark` represents a marks path that has to be 
  !! followed. 
  !! Otherwise, `mark` is a unique mark that has to be searched throughout the file.  
  module procedure xread_int
  module procedure xread_real64
  module procedure xread_char
  module procedure xread_int_alloc
  module procedure xread_real64_alloc
  module procedure xread_char_alloc
end interface

contains

!-----------------------------------------------------------------------
! xread_int (xread)
!-----------------------------------------------------------------------
function xread_int(id, mark, d, rew) result(res)
!! Read an integer scalar.
integer,           intent(in)  :: id   !! File unit
character(*),      intent(in)  :: mark !! Mark path to be followed or unique mark to be searched.
integer,           intent(out) :: d    !! Output variable.
logical, optional, intent(in)  :: rew  !! Whether to rewind the file after read.
integer :: res
character(maxpath) :: str

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
backspace(id)
read (id, *, iostat=res) str, d
if (res /= 0) call error('(xread_int/read) error while reading mark "'//&
trim(mark)//'", #'//trim(string(res)))
call info('read <'//trim(adjustlt(mark))//'>: '//string(d))
end function

!-----------------------------------------------------------------------
! xread_real64 (xread)
!-----------------------------------------------------------------------
function xread_real64(id, mark, d, rew) result(res)
!! Read a real64 scalar.
integer,           intent(in)  :: id   !! File unit
character(*),      intent(in)  :: mark !! Mark path to be followed or unique mark to be searched.
real(real64),      intent(out) :: d    !! Output variable.
logical, optional, intent(in)  :: rew  !! Whether to rewind the file after read.
integer :: res
character(maxpath) :: str

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
backspace(id)
read (id, *, iostat=res) str, d
if (res /= 0) call error('(xread_real64/read) error while reading mark "'//&
trim(mark)//'", #'//trim(string(res)))
call info('read <'//trim(adjustlt(mark))//'>: '//string(d))
end function

!-----------------------------------------------------------------------
! xread_char (xread)
!-----------------------------------------------------------------------
function xread_char(id, mark, d, rew) result(res)
!! Read a character scalar.
integer,           intent(in)  :: id   !! File unit
character(*),      intent(in)  :: mark !! Mark path to be followed or unique mark to be searched.
character(*),      intent(out) :: d    !! Output variable.
logical, optional, intent(in)  :: rew  !! Whether to rewind the file after read.
integer :: res

if (search_path(id, mark, rew) /= 0) then; res = -1; return; end if
read (id, '(a)', iostat=res) d
if (res /= 0) call error('(xread_char/read) error while reading mark "'//&
trim(mark)//'", #'//trim(string(res)))
call info('read <'//trim(adjustlt(mark))//'>: '//'"'//trim(d)//'"')
end function

!-----------------------------------------------------------------------
! xread_int_alloc (xread)
!-----------------------------------------------------------------------
function xread_int_alloc(id, mark, d, rew) result(res)
!! Read an integer allocatable array.
integer,              intent(in)  :: id   !! File unit
character(*),         intent(in)  :: mark !! Mark path to be followed or unique mark to be searched.
integer, allocatable, intent(out) :: d(:) !! Output variable.
logical, optional,    intent(in)  :: rew  !! Whether to rewind the file after read.
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
! xread_real64_alloc (xread)
!-----------------------------------------------------------------------
function xread_real64_alloc(id, mark, d, rew) result(res)
!! Read a real64 allocatable array.
integer,                   intent(in)  :: id   !! File unit
character(*),              intent(in)  :: mark !! Mark path to be followed or unique mark to be searched.
real(real64), allocatable, intent(out) :: d(:) !! Output variable.
logical, optional,         intent(in)  :: rew  !! Whether to rewind the file after read.
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
! xread_char_alloc (xread)
!-----------------------------------------------------------------------
function xread_char_alloc(id, mark, d, rew) result(res)
!! Read a character allocatable array.
integer,                   intent(in)  :: id   !! File unit
character(*),              intent(in)  :: mark !! Mark path to be followed or unique mark to be searched.
character(*), allocatable, intent(out) :: d(:) !! Output variable.
logical, optional,         intent(in)  :: rew  !! Whether to rewind the file after read.
integer :: res
character(maxpath) :: str, cad
character(4096) :: longstr
integer :: n, i

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
! xlist
!-----------------------------------------------------------------------
subroutine xlist(id, path, val, rew)
!! List all the marks contained in a mark path.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: maxpath, xlist`  
!! `implicit none`  
!! `character(maxpath), allocatable :: val(:)`  
!! `open(10, file='data.xml')`  
!! `call xlist(10, '/data', val, rew=.true.)`  
!! `close(10)`  
!! `end program`  
integer,                   intent(in)  :: id     !! File unit
character(*),              intent(in)  :: path   !! Mark path to be followed.
character(*), allocatable, intent(out) :: val(:) !! Output variable.
logical, optional,         intent(in)  :: rew    !! Whether to rewind the file after read.
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
! is_alphanumeric
!-----------------------------------------------------------------------
function is_alphanumeric(c) result(res)
!! Checks whether the first character is alphanumeric (alphabetic or numeric character).  
!! @note `iachar` only considers the first character of `c`.
character(*) :: c   !! Character to test.
logical      :: res !! `.true.` when `c` is alphanumeric.

if ((iachar('a') <= iachar(c) .and. iachar(c) <= iachar('z')) .or. &
    (iachar('A') <= iachar(c) .and. iachar(c) <= iachar('Z')) .or. &
    (iachar('0') <= iachar(c) .and. iachar(c) <= iachar('9'))) then
  res = .true.
else
  res = .false.
end if
end function

!--------------------------------------------------------------------
! aword
!--------------------------------------------------------------------
function aword(cad) result(res)
!! Gets the first word of a string removing the trailing non-alphanumeric characters.  
!! @note Please use `trim(aword(cad))`.
character(*), intent(in) :: cad !! Input string.
character(maxpath)       :: res !! First word in the string.
integer :: ios

read(cad, *, iostat=ios) res
if (ios /= 0) call error('aword/read, #'//trim(string(ios)))
res = adjustl(res)
do while (.not. is_alphanumeric(res(len_trim(res):len_trim(res))))
  res = res(1:len_trim(res)-1)
end do
end function

!-----------------------------------------------------------------------
! search_path
!-----------------------------------------------------------------------
function search_path(id, mark, rew) result(res)
!! Searchs a mark or a path; if mark(1:1) is alphanumeric, 
!! a mark is searched with procedure search_mark; otherwise, 'mark' is 
!! considered a path of marks separated by mark(1:1).  
!!
!! @note `0` if the mark is found; `iostat_end` if end-of-file is found; `non-zero` otherwise.
integer,      intent(in)      :: id
character(*), intent(in)      :: mark
logical, optional, intent(in) :: rew
integer                       :: res
character(maxpath) :: part, path
character(1) :: sep
integer :: p

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
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! search_once
!-----------------------------------------------------------------------
function search_once(id, mark) result(res)
!! Searches for a mark only once.  
!!
!! @note  0 if the mark is found; iostat_end if end-of-file is found; non-zero otherwise
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
! search_mark
!-----------------------------------------------------------------------
function search_mark(id, mark) result(res)
!! Searches for `<mark` in the unit id.  
!!
!! @note 0 if the mark is found; iostat_end if end-of-file is found; non-zero otherwise
integer,      intent(in) :: id
character(*), intent(in) :: mark
integer                  :: res

do
  res = search_once(id, '<'//trim(adjustlt(lcase(mark))))
  if (res == iostat_end .or. res == 0) return !end-of-file or mark found
end do
end function

!-----------------------------------------------------------------------
! search_close_mark
!-----------------------------------------------------------------------
function search_close_mark(id, mark) result(res)
  !! Searchs a close mark.  
integer,      intent(in) :: id
character(*), intent(in) :: mark

integer                  :: res

character(maxpath) :: str
integer :: n

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
! last_mark
!-----------------------------------------------------------------------
function last_mark(path) result(mark)
  !! Finds the last mark of a path.  
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

end module
