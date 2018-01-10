module module_convers
!-----------------------------------------------------------------------
! Module for type conversion
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 18/01/2012
!
! PUBLIC CONSTANTS:
! default_format_*: recommended data descriptors for default output
!
! PUBLIC PROCEDURES:
!   string: converts a numeric or logical variable (scalar or array) to character
!   int: extends int to convert a character to int
!   int_alloc: extends int to convert a character to an integer allocatable array
!   real: extends real to convert a character to real
!   dble: extends dble to convert a character to double_precision
!   dble_alloc: extends dble to convert a character to a real64 allocatable array
!   word_count: counts the number of words in a string
!   int_count: counts the number of integer constants in a string
!   dble_count: counts the number of real64 constants in a string
!   lcase: converts a string to lower case
!   adjustlt: extension of adjustl for removing left tab characters (char(9))
!   replace: subroutine to replace in a string a substring by another one
!   freplace: function that returns a string replacing a substring by another one
!   word: gets the first word of a string
!   trim: returns the begining of a string before the first/last appearance of separator
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error
implicit none

!Constants
character(*), parameter :: default_format_integer = 'I12'    ! default integer format
character(*), parameter :: default_format_real    = 'E15.7'  ! default real format
character(*), parameter :: default_format_real64  = 'E25.16' ! default real64 format

!Private procedures
private :: string_int, string_real, string_dbl, string_log, string_char, int_char, &
           int_char_alloc, real_char, dble_char, string_int_array, string_dble_array, &
           trim_prv

!Interfaces
interface string;  module procedure string_int;        end interface
interface string;  module procedure string_real;       end interface
interface string;  module procedure string_dbl;        end interface
interface string;  module procedure string_log;        end interface
interface string;  module procedure string_char;       end interface
interface string;  module procedure string_int_array;  end interface
interface string;  module procedure string_dble_array; end interface
interface int;     module procedure int_char;          end interface
interface int;     module procedure int_char_alloc;    end interface
interface real;    module procedure real_char;         end interface
interface dble;    module procedure dble_char;         end interface
interface trim;    module procedure trim_prv;          end interface

contains

!--------------------------------------------------------------------
! string: converts a numeric variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_int(i) result(res)
integer, intent(in) :: i
character(maxpath) :: res, iom, aux
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) i
if (ios /= 0) then
  write(unit=aux, fmt=*) ios
  call error('string_int/write, #'//trim(aux)//': '//trim(iom))
endif
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string: converts a numeric variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_real(x) result(res)
real, intent(in) :: x
character(maxpath) :: res, iom
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) x
if (ios /= 0) call error('string_real/write, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string: converts a numeric variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_dbl(d) result(res)
real(real64), intent(in) :: d
character(maxpath) :: res, iom
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) d
if (ios /= 0) call error('string_dbl/write, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string: converts a logical variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_log(l) result(res)
logical, intent(in) :: l
character(maxpath) :: res, iom
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) l
if (ios /= 0) call error('string_log/write, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string: converts a char variable to character (erase quotation)
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_char(l) result(res)
character(*), intent(in) :: l
character(maxpath) :: res, iom
integer :: ios

read(unit=l, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('string_char/read, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string: converts a integer array variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_int_array(x) result(res)
integer, intent(in) :: x(:)
character(13*size(x,1)) :: res
integer :: i

res = ' '
do i = 1, size(x,1)
  res = trim(res)//' '//trim(adjustlt(string(x(i))))
end do  
end function

!--------------------------------------------------------------------
! string: converts a real64 array variable to character
! NOTE: please use trim(string(num))
!--------------------------------------------------------------------
function string_dble_array(x) result(res)
real(real64), intent(in) :: x(:)
character(26*size(x,1)) :: res
integer :: i

res = ' '
do i = 1, size(x,1)
  res = trim(res)//' '//trim(adjustlt(string(x(i))))
end do  
end function

!-----------------------------------------------------------------------
! int: extends int to convert a character to int
!-----------------------------------------------------------------------
function int_char(cad) result(res)
character(*), intent(in) :: cad
integer                      :: res
character(maxpath) :: iom
integer :: ios

read(unit=cad, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('int_char/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! int: extends int to convert a character to int
!-----------------------------------------------------------------------
function int_char_alloc(cad) result(res)
character(*), dimension(:), allocatable, intent(in) :: cad
integer, dimension(size(cad,1)) :: res
character(maxpath) :: iom
integer :: ios, i

do i = 1, size(cad,1)
  read(unit=cad(i), fmt=*, iostat=ios, iomsg=iom) res(i)
enddo
if (ios /= 0) call error('int_char_alloc/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! int_alloc: extends int to convert a character to an integer allocatable array
!-----------------------------------------------------------------------
function int_alloc(str) result(res)
character(*), intent(in) :: str
integer, allocatable  :: res(:)
character(len(str)) :: tmp
character(maxpath) :: cad
integer :: ios, i, p, n

!allocation
n = word_count(str)
allocate(res(n), stat = ios, errmsg = cad)
if (ios /= 0) call error('(module_convers/int_alloc) unable to allocate variable: '//trim(cad))
!reading
tmp = adjustlt(str)
do i = 1, n
  p = index(tmp, ' ')
  read(tmp(1:p-1), *, iostat = ios) res(i)
  if (ios /= 0) call error('(module_convers/int_alloc) unable to read '//trim(adjustl(string(i)))//'th data: '//trim(tmp(p+1:)))
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
end function

!-----------------------------------------------------------------------
! real: extends real to convert a character to real
!-----------------------------------------------------------------------
function real_char(cad) result(res)
character(*), intent(in) :: cad
real :: res
character(maxpath) :: iom
integer :: ios

read(unit=cad, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('real_char/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! dble: extends dble to convert a character to double precision
!-----------------------------------------------------------------------
function dble_char(cad) result(res)
character(*), intent(in) :: cad
real(real64) :: res
character(maxpath) :: iom
integer :: ios

read(unit=cad, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('dble/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! dble_alloc: extends dble to convert a character to a real64 allocatable array
!-----------------------------------------------------------------------
function dble_alloc(str) result(res)
character(*), intent(in) :: str
real(real64), allocatable  :: res(:)
character(len(str)) :: tmp
character(maxpath) :: cad
integer :: ios, i, p, n

!allocation
n = word_count(str)
allocate(res(n), stat = ios, errmsg = cad)
if (ios /= 0) call error('(module_convers/int_alloc) unable to allocate variable: '//trim(cad))
!reading
tmp = adjustlt(str)
do i = 1, n
  p = index(tmp, ' ')
  read(tmp(1:p-1), *, iostat = ios) res(i)
  if (ios /= 0) call error('(module_convers/int_alloc) unable to read '//trim(adjustl(string(i)))//'th data: '//trim(tmp(p+1:)))
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
end function

!-----------------------------------------------------------------------
! word_count: counts the number of words in a string separated by sep
!-----------------------------------------------------------------------
function word_count(str, sep) result(res)
character(*),           intent(in) :: str
character(*), optional, intent(in) :: sep
character(len(str))                :: tmp
integer :: res, p

res = 0
tmp = adjustlt(str)
do 
  if (len_trim(tmp) == 0) return !empty string
  if(present(sep)) then
    p = index(trim(tmp), sep)  
  else
    p = index(trim(tmp), ' ')  
  endif
  if (p == 0) then !separator not found
    if (len_trim(tmp) > 0) res = res + 1
    return
  end if
  if (len_trim(tmp(1:p-1)) == 0) then !there is nothing before separator
    tmp = adjustlt(tmp(p+1:len_trim(tmp)))
    cycle
  end if
  res = res + 1 
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
  
end function

!-----------------------------------------------------------------------
! int_count: counts the number of integers in a string
!-----------------------------------------------------------------------
function int_count(str) result(res)
character(*), intent(in) :: str
character(len(str)) :: tmp
integer :: res, p, ios
integer :: x

res = 0
tmp = adjustlt(str)
do 
  p = index(tmp, ' ')
  if (tmp(1:p-1) == ' ') return 
  read(tmp(1:p-1),*, iostat = ios) x
  if (ios == 0) res = res + 1 
  tmp = trim(adjustlt(tmp(p+1:len_trim(tmp))))
end do
  
end function

!-----------------------------------------------------------------------
! dble_count: counts the number of numbers in a string
!-----------------------------------------------------------------------
function dble_count(str) result(res)
character(*), intent(in) :: str
character(len(str)) :: tmp
integer :: res, ios, p
real(real64) :: x

res = 0
tmp = adjustlt(str)
do 
  p = index(tmp, ' ')
  if (tmp(1:p-1) == ' ') return 
  read(tmp(1:p-1),*, iostat = ios) x
  if (ios == 0) res = res + 1 
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
  
end function

!-----------------------------------------------------------------------
! lcase: converts a string to lower case
!-----------------------------------------------------------------------
function lcase(str) result(res)
character(*), intent(in) :: str
character(len(str)) :: res
integer :: diff, i

res = str
diff = ichar('A') - ichar('a')
do i = 1, len_trim(str)
  if (str(i:i) < 'A' .or. str(i:i) > 'Z') cycle
  res(i:i) = char(ichar(str(i:i)) - diff)
end do
end function

!-----------------------------------------------------------------------
! adjustlt: extension of adjustl to remove tab characters (char(9))
!-----------------------------------------------------------------------
function adjustlt(string) result(res)
character(*), intent(in) :: string
character(len(string)) :: res

res = string
if (len_trim(string) <= 0) return
do while ((res(1:1) == char(9) .or. res(1:1) == ' ') .and. len_trim(res)>0)
  res = res(2:len_trim(res))
enddo
end function

!-----------------------------------------------------------------------
! replace: replace in formula cad1 by cad2
!-----------------------------------------------------------------------
subroutine replace(formula, cad1, cad2)
character(*), intent(inout) :: formula
character(*), intent(in) :: cad1, cad2
integer :: pos

pos = index(trim(formula), cad1)
do while (pos>0)
  formula = formula(1:pos-1)//cad2//formula(pos+len(cad1):len_trim(formula))
  pos = index(trim(formula), cad1)
enddo
end subroutine

!-----------------------------------------------------------------------
! freplace: replace in formula cad1 by cad2
!-----------------------------------------------------------------------
function freplace(formula, cad1, cad2) result(res)
character(*), intent(in) :: formula, cad1, cad2
character(len(formula)) :: res

res = formula
call replace(res, cad1,cad2)
end function

!--------------------------------------------------------------------
! word: gets the n-th word of a string
! NOTE: please use trim(word(cad))
!--------------------------------------------------------------------
function word(str, pos) result(res)
character(*),      intent(in) :: str
integer, optional, intent(in) :: pos
character(len=len(str)) :: res, tmp
integer :: ios, n, p, i

n = 1
if (present(pos)) n = pos
if (n > word_count(str)) call error('(module_convers/word) requested position is beyong string length')
!search n-th word
tmp = adjustlt(str)
do i = 1, n-1
  p = index(tmp, ' ')
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
!rad n-th word
read(tmp, *, iostat=ios) res
if (ios /= 0) call error('(module_convers/word) unable to read from string, #'//trim(string(ios)))
res = adjustlt(res)
end function

!-----------------------------------------------------------------------
! trim: returns the beginning of a string before the first/last appearance of separator
!-----------------------------------------------------------------------
function trim_prv(str, sep, back) result(res)
character(*),      intent(in) :: str
character(*),      intent(in) :: sep
logical, optional, intent(in) :: back
character(len(str)) :: res
integer :: pos

pos = index(str, sep, back = back)
if (pos == 0) then
  res = str
else
  res = str(1:pos-1)
endif
end function

end module
