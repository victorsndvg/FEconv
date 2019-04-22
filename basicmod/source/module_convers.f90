module module_convers_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for type conversion.
!
! PUBLIC CONSTANTS:
!   default_format_integer: Default format for integer variables
!   default_format_real: Default format for real variables
!   default_format_real64: Default format for real64 format
!
! PUBLIC PROCEDURES:
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
use module_compiler_dependant_bmod, only: real64
use module_os_dependant_bmod, only: maxpath
use module_report_bmod, only: error
use module_alloc_common_bmod, only: csize
implicit none

 !Constants
character(*), parameter :: default_format_integer = 'I12'    !! Default format for integer variables.
character(*), parameter :: default_format_real    = 'E15.7'  !! Default format for real variables.
character(*), parameter :: default_format_real64  = 'E25.16' !! Default format for real64 format.

 !Private procedures
private :: string_int, string_real, string_dbl, string_log, string_char, int_char, &
           int_char_vec, real_char, dble_char, string_int_vec, string_dble_vec, &
           trim_prv, is_int_sca_prv, is_int_vec_prv, is_real_sca_prv, is_real_vec_prv, &
           is_dble_sca_prv, is_dble_vec_prv

interface string
  !! Converts the input variable to character. If the input variable is already character,
  !! returns its value without quotation marks.  
  !!
  !! @note It is recommended to use `trim(string(argument))` to avoid trailing blanks.
  module procedure string_int
  module procedure string_real
  module procedure string_dbl
  module procedure string_log
  module procedure string_char
  module procedure string_int_vec
  module procedure string_dble_vec
end interface

interface int
  !! Converts a character scalar or array into integer. If the input variable is an array, the result is an array of the same
  !! shape.  
  !!
  !! @warning If the input variable cannot be converted to integer, an error is raised.
  module procedure int_char
  module procedure int_char_vec
end interface

interface real
  !! Converts a character scalar to real.
  !!
  !! @warning If the input variable cannot be converted to real, an error is raised.
  module procedure real_char
end interface

interface dble
  !! Converts a character scalar to real64.
  !!
  !! @warning If the input variable cannot be converted to real64, an error is raised.
  module procedure dble_char
end interface

interface is_int
  !! Checks whether a character scalar or array contains integer(s).
  !!
  !! @note The ideas for this procedure were taken from
  !! [this Intel forum](https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/285098).  
  !!
  !! @note The scalar version has not been declared _elemental_ for compatibility with old compilers.
  !! Thus, the array version had to be implemented.  
  module procedure is_int_sca_prv
  module procedure is_int_vec_prv
end interface

interface is_real
  !! Checks whether a character scalar or array contains real(s).
  !!
  !! @note The ideas for this procedure were taken from
  !! [this Intel forum](https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/285098).  
  !!
  !! @note The scalar version has not been declared _elemental_ for compatibility with old compilers.
  !! Thus, the array version had to be implemented.  
  module procedure is_real_sca_prv
  module procedure is_real_vec_prv
end interface

interface is_dble
  !! Checks whether a character scalar or array contains double(s).
  !!
  !! @note The ideas for this procedure were taken from
  !! [this Intel forum](https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/285098).  
  !!
  !! @note The scalar version has not been declared _elemental_ for compatibility with old compilers.
  !! Thus, the array version had to be implemented.  
  module procedure is_dble_sca_prv
  module procedure is_dble_vec_prv
end interface

interface trim
  !! Returns the begining of a string before the first/last appearance of a separator.
  module procedure trim_prv
end interface

contains

!--------------------------------------------------------------------
! string_int (string)
!--------------------------------------------------------------------
function string_int(i) result(res)
!! Converts an integer scalar to character.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `print*, trim(string(5))`  
!! `end program`  
integer, intent(in) :: i  !! Variable to convert to character.
character(maxpath) :: res !! Value of the input variable as character.
character(maxpath) :: iom, aux
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) i
if (ios /= 0) then
  write(unit=aux, fmt=*) ios
  call error('module_convers/string_int/write, #'//trim(aux)//': '//trim(iom))
endif
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string_real (string)
!--------------------------------------------------------------------
function string_real(x) result(res)
!! Converts a real scalar to character.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `print*, trim(string(5.56))`  
!! `end program`  
real, intent(in) :: x     !! Variable to convert to character.
character(maxpath) :: res !! Value of the input variable as character.
character(maxpath) :: iom
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) x
if (ios /= 0) call error('module_convers/string_real/write, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string_dbl (string)
!--------------------------------------------------------------------
function string_dbl(d) result(res)
!! Converts a real64 scalar to character.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `print*, trim(string(5.56_real64))`  
!! `end program`  
real(real64), intent(in) :: d !! Variable to convert to character.
character(maxpath) :: res     !! Value of the input variable as character.
character(maxpath) :: iom
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) d
if (ios /= 0) call error('module_convers/string_dbl/write, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string_log (string)
!--------------------------------------------------------------------
function string_log(l) result(res)
!! Converts a logical scalar to character.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `print*, trim(string(.true.))`  
!! `end program`  
logical, intent(in) :: l  !! Variable to convert to character.
character(maxpath) :: res !! Value of the input variable as character.
character(maxpath) :: iom
integer :: ios

write(unit=res, fmt=*, iostat=ios, iomsg=iom) l
if (ios /= 0) call error('module_convers/string_log/write, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string_char (string)
!--------------------------------------------------------------------
function string_char(c) result(res)
!! Converts a character scalar to character (erase quotation).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `character(20) :: c = '''testing.file.txt'''`  
!! `print*, trim(string(c))`  
!! `end program`  
character(*), intent(in) :: c !! Variable to convert to character.
character(maxpath) :: res     !! Contains the value of the input variable without quotation marks.
character(maxpath) :: iom
integer :: ios

read(unit=c, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('module_convers/string_char/read, #'//trim(string(ios))//': '//trim(iom))
res = adjustl(res)
end function

!--------------------------------------------------------------------
! string_int_vec (string)
!--------------------------------------------------------------------
function string_int_vec(x) result(res)
!! Converts a integer array to character scalar with blank-separated words.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `integer :: x(5) = [1,-3,0,100,4]`  
!! `print*, trim(string(x))`  
!! `end program`  
integer, intent(in)     :: x(:) !! Integer array to convert.
character(13*size(x,1)) :: res  !! Value of the input variable as character; the length of each integer is assumed to be <= 12
                                !! (see variable `default_format_integer` in the `module_convers` documentation.
integer :: i

res = ' '
do i = 1, size(x,1)
  res = trim(res)//' '//trim(adjustlt(string(x(i))))
end do
res = adjustlt(res)
end function

!--------------------------------------------------------------------
! string_dbl_vec (string)
!--------------------------------------------------------------------
function string_dble_vec(x) result(res)
!! Converts a real64 array to character scalar with blank-separated words.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: string`  
!! `implicit none`  
!! `real(real64) :: x(5) = [1.,-3.,0.,100.,4.]`  
!! `print*, trim(string(x))`  
!! `end program`  
real(real64), intent(in) :: x(:) !! Real64 array to convert.
character(26*size(x,1))  :: res  !! Value of the real64 variable as character; the length of each real64 is assumed to be <= 25
                                 !! (see variable `default_format_real64` in the `module_convers` documentation.
integer :: i

res = ' '
do i = 1, size(x,1)
  res = trim(res)//' '//trim(adjustlt(string(x(i))))
end do
res = adjustlt(res)
end function

!-----------------------------------------------------------------------
! int_char (int)
!-----------------------------------------------------------------------
function int_char(cad) result(res)
!! Converts a character scalar to integer.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: int`  
!! `implicit none`  
!! `print*, int('5')+2`  
!! `end program`  
character(*), intent(in) :: cad  !! Character variable to be converted.
integer                  :: res  !! Integer value, if the conversion was possible.
character(maxpath) :: iom
integer :: ios

read(unit=cad, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('module_convers/int_char/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! int_char_vec (int)
!-----------------------------------------------------------------------
function int_char_vec(cad) result(res)
!! Converts a character array to an integer array of the same shape.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: int`  
!! `implicit none`  
!! `character(12) :: cad(3) = ['1','234','-99']`  
!! `integer :: res(3)`  
!! `res = int(cad)`  
!! `print*, res`  
!! `end program`  
!!
!! @note An error is raised if a component cannot converted to integer or the array has no components.
character(*), dimension(:), intent(in) :: cad !! Character array to be converted.
integer,      dimension(size(cad,1))   :: res !! Integer array of the same shape than the input one.
integer :: i

if (csize(char1=cad, d=1) <= 0) call error('(module_convers/int_char_vec), input array seems to have no componentes.')
do i = 1, size(cad,1)
  if (.not. is_int(cad(i))) call error('(module_convers/int_char_vec), a component cannot be converted to integer, '//&
  &trim(string(cad(i))))
  res(i) = int(cad(i))
enddo
end function

!-----------------------------------------------------------------------
! int_alloc
!-----------------------------------------------------------------------
subroutine int_alloc(str, res)
!! Extract integers from a character scalar with blank-separated words to create an integer array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: int_alloc`  
!! `implicit none`  
!! `character(12) :: cad = '1 234 -99 numbers'`  
!! `integer, allocatable :: res(:)`  
!! `call int_alloc(cad, res)`  
!! `end program`  
character(*),         intent(in)  :: str    !! Character array to be converted.
integer, allocatable, intent(out) :: res(:) !! Integer array with the converted integers.
character(len(str)) :: tmp
character(maxpath) :: cad
integer :: ios, i, p, n

!allocation
n = int_count(str)
allocate(res(n), stat = ios, errmsg = cad)
if (ios /= 0) call error('(module_convers/int_alloc) unable to allocate variable: '//trim(cad))
!reading
tmp = adjustlt(str)
i = 1
do while (i <= n)
  p = index(tmp, ' ')
  if (is_int(tmp(1:p-1))) then
    res(i) = int(tmp(1:p-1))
    i = i+1
  end if
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
end subroutine

!-----------------------------------------------------------------------
! real_char (real)
!-----------------------------------------------------------------------
function real_char(cad) result(res)
!! Converts a scalar character to real.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real`  
!! `implicit none`  
!! `print*, real('5.4')+2.`  
!! `end program`  
character(*), intent(in) :: cad !! Character variable to be converted.
real                     :: res !! Real value, if the conversion was possible.
character(maxpath) :: iom
integer :: ios

read(unit=cad, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('module_convers/real_char/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! dble_char (dble)
!-----------------------------------------------------------------------
function dble_char(cad) result(res)
!! Converts a character scalar to real64.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dble`  
!! `implicit none`  
!! `print*, dble('5.4_real64')+2.`  
!! `end program`  
character(*), intent(in) :: cad !! Character variable to be converted.
real(real64)             :: res !! Real64 value, if the conversion was possible.
character(maxpath) :: iom
integer :: ios

read(unit=cad, fmt=*, iostat=ios, iomsg=iom) res
if (ios /= 0) call error('module_convers/dble/read, #'//trim(string(ios))//': '//trim(iom))
end function

!-----------------------------------------------------------------------
! dble_alloc
!-----------------------------------------------------------------------
subroutine dble_alloc(str, res)
!! Extract doubles from a character scalar with blank-separated words to create an real64 array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dble_alloc`  
!! `implicit none`  
!! `character(12) :: cad = '1 234 -99 numbers'`  
!! `real(real64), allocatable :: res(:)`  
!! `call dble_alloc(cad, res)`  
!! `end program`  
character(*),              intent(in)  :: str    !! Character array to be converted.
real(real64), allocatable, intent(out) :: res(:) !! Real64 array with the converted doubles.
character(len(str)) :: tmp
character(maxpath) :: cad
integer :: ios, i, p, n

!allocation
n = dble_count(str)
allocate(res(n), stat = ios, errmsg = cad)
if (ios /= 0) call error('(module_convers/dble_alloc) unable to allocate variable: '//trim(cad))
!reading
tmp = adjustlt(str)
do while (i <= n)
  p = index(tmp, ' ')
  if (is_dble(tmp(1:p-1))) then
    res(i) = dble(tmp(1:p-1))
    i = i+1
  end if
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
end subroutine

!-----------------------------------------------------------------------
! is_int_sca_prv (is_int)
!-----------------------------------------------------------------------
function is_int_sca_prv(str) result(res)
!! Checks whether a character scalar contains an integer.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_int`  
!! `implicit none`  
!! `character(120) :: cad = '1 1' !not an integer`  
!! `print*, is_int(cad)`  
!! `end program`  
character(*), intent(in) :: str !! Character scalar to be analyzed.
logical                  :: res !! `.true.` if the string is analyzed as integer; `.false.` otherwise.
integer :: i, ios

read(str, *, iostat = ios) i
res = (ios == 0) .and. (verify(trim(adjustlt(str)), '+-0123456789') == 0)
end function

!-----------------------------------------------------------------------
! is_int_vec_prv (is_int)
!-----------------------------------------------------------------------
function is_int_vec_prv(str) result(res)
!! Checks whether a character array contains integers.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_int`  
!! `implicit none`  
!! `character(120) :: cad(2) = ['1 1', '123']`  
!! `print*, is_int(cad)`  
!! `end program`  
character(*), intent(in) :: str(:) !! Character array to be analyzed.
logical                  :: res !! `.true.` if every component is analyzed as integer; `.false.` otherwise.
integer :: i

res = .false.
if (csize(char1=str, d=1) <= 0) return
res = all([(is_int(str(i)), i=1, csize(char1=str, d=1))])
end function

!-----------------------------------------------------------------------
! is_real_sca_prv (is_real)
!-----------------------------------------------------------------------
function is_real_sca_prv(str) result(res)
!! Checks whether a character scalar contains a real.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_real`  
!! `implicit none`  
!! `character(120) :: cad = '1.1e1'`  
!! `print*, is_real(cad)`  
!! `end program`  
character(*), intent(in) :: str !! Character scalar to be analyzed.
logical                  :: res !! `.true.` if the string is analyzed as real; `.false.` otherwise.
integer :: ios
real :: x

read(str, *, iostat = ios) x
res = (ios == 0) .and. (verify(trim(adjustlt(str)), '+-.0123456789eE') == 0)
end function

!-----------------------------------------------------------------------
! is_real_vec_prv (is_real)
!-----------------------------------------------------------------------
function is_real_vec_prv(str) result(res)
!! Checks whether a character array contains reals.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_real`  
!! `implicit none`  
!! `character(120) :: cad(2) = ['1 1 ', '1.23']`  
!! `print*, is_real(cad)`  
!! `end program`  
character(*), intent(in) :: str(:) !! Character array to be analyzed.
logical                  :: res !! `.true.` if every component is analyzed as real; `.false.` otherwise.
integer :: i

res = .false.
if (csize(char1=str, d=1) <= 0) return
res = all([(is_real(str(i)), i=1, csize(char1=str, d=1))])
end function

!-----------------------------------------------------------------------
! is_dble_sca_prv (is_dble)
!-----------------------------------------------------------------------
function is_dble_sca_prv(str) result(res)
!! Checks whether a character scalar contains a double.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_dble`  
!! `implicit none`  
!! `character(120) :: cad = '1.1D1'`  
!! `print*, is_dble(cad)`  
!! `end program`  
character(*), intent(in) :: str !! Character scalar to be analyzed.
logical                  :: res !! `.true.` if the string is analyzed as double; `.false.` otherwise.
integer :: ios
real(real64) :: x

read(str, *, iostat = ios) x
res = (ios == 0) .and. (verify(trim(adjustlt(str)), '+-.0123456789eEdD') == 0)
end function

!-----------------------------------------------------------------------
! is_dble_vec_prv (is_dble)
!-----------------------------------------------------------------------
function is_dble_vec_prv(str) result(res)
!! Checks whether a character array contains doubles.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: is_dble`  
!! `implicit none`  
!! `character(120) :: cad(2) = ['1 1 ', '1.D3']`  
!! `print*, is_dble(cad)`  
!! `end program`  
character(*), intent(in) :: str(:) !! Character array to be analyzed.
logical                  :: res !! `.true.` if every component is analyzed as double; `.false.` otherwise.
integer :: i

res = .false.
if (csize(char1=str, d=1) <= 0) return
res = all([(is_dble(str(i)), i=1, csize(char1=str, d=1))])
end function

!-----------------------------------------------------------------------
! word_count
!-----------------------------------------------------------------------
function word_count(str, sep) result(res)
!! Counts the number of words in a string separated by a separator.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: word_count`  
!! `implicit none`  
!! `character(120) :: cad = 'a number inside'`  
!! `print*, word_count(cad,' ')`  
!! `end program`  
character(*),           intent(in) :: str !! String to be analyzed.
character(*), optional, intent(in) :: sep !! Separator.
integer                            :: res !! Number of words in a string separated by a separator.
character(len(str)) :: tmp, wrd, part
integer :: p, lsep, lquote, ios

res = 0
tmp = adjustlt(str)
if (present(sep)) then
  lsep = len_trim(sep)
  do
    if (len_trim(tmp) == 0) return !empty string
    p = index(trim(tmp), sep)
    if (p == 0) then !separator not found
      if (len_trim(tmp) > 0) res = res + 1
      return
    end if
    if (len_trim(tmp(1:p-1)) == 0) then !there is nothing before separator
      tmp = adjustlt(tmp(p+lsep:len_trim(tmp)))
      cycle
    end if
    res = res + 1
    tmp = adjustlt(tmp(p+lsep:len_trim(tmp)))
  end do
else
  do
    if (len_trim(tmp) == 0) return !empty string
    if (tmp(1:1) == '"' .or. tmp(1:1) == "'") then
      lquote = 2
    else
      lquote = 0
    end if
    read(tmp, *, iostat=ios) wrd
    if (ios /= 0) call error('(module_convers/word_count) unable to read from string '//trim(tmp)//', #'//trim(string(ios)))
    if (len_trim(wrd) == 0) return !empty string
    wrd = adjustlt(wrd)
    tmp = adjustlt(tmp)
    tmp = adjustlt(tmp(len_trim(wrd)+lquote+1:))
    !detect slash
    do while (tmp(1:1) == '/')
      tmp = tmp(2:)
      read(tmp, *, iostat=ios) part
      if (ios /= 0) call error('(module_convers/word_count) unable to read from string '//trim(tmp)//', #'//trim(string(ios)))
      part = adjustlt(part)
      tmp  = adjustlt(tmp)
      tmp  = adjustlt(tmp(len_trim(part)+1:))
    end do
    !detect comma(s)
    do while (tmp(1:1) == ',')
      tmp = adjustlt(tmp(2:))
    end do
    res = res + 1
  end do
end if
end function

!-----------------------------------------------------------------------
! int_count
!-----------------------------------------------------------------------
function int_count(str) result(res)
!! Counts in a string the number of blank-separated words that can be converted to integer.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: int_count`  
!! `implicit none`  
!! `character(120) :: cad = 'just 1 integer'`  
!! `print*, int_count(cad)`  
!! `end program`  
character(*), intent(in) :: str !! String to be analyzed.
integer                  :: res !! Number of blank-separated words that can be converted to integer in the string.
character(len(str)) :: tmp
integer :: p, ios

res = 0
tmp = adjustlt(str)
do
  p = index(tmp, ' ')
  if (tmp(1:p-1) == ' ') return
  if (is_int(tmp(1:p-1))) res = res + 1
  tmp = trim(adjustlt(tmp(p+1:len_trim(tmp))))
end do
end function

!-----------------------------------------------------------------------
! real_count
!-----------------------------------------------------------------------
function real_count(str) result(res)
!! Counts in a string the number of blank-separated words that can be converted to real64.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real_count`  
!! `implicit none`  
!! `character(120) :: cad = 'just 1. real'`  
!! `print*, real_count(cad)`  
!! `end program`  
character(*), intent(in) :: str !! String to be analyzed.
integer                  :: res !! Number of blank-separated words that can be converted to real64 in the string.
character(len(str)) :: tmp
integer :: ios, p

res = 0
tmp = adjustlt(str)
do
  p = index(tmp, ' ')
  if (tmp(1:p-1) == ' ') return
  if (is_real(tmp(1:p-1))) res = res + 1
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
end function

!-----------------------------------------------------------------------
! dble_count
!-----------------------------------------------------------------------
function dble_count(str) result(res)
!! Counts in a string the number of blank-separated words that can be converted to real64.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dble_count`  
!! `implicit none`  
!! `character(120) :: cad = 'just 1. double'`  
!! `print*, dble_count(cad)`  
!! `end program`  
character(*), intent(in) :: str !! String to be analyzed.
integer                  :: res !! Number of blank-separated words that can be converted to real64 in the string.
character(len(str)) :: tmp
integer :: ios, p

res = 0
tmp = adjustlt(str)
do
  p = index(tmp, ' ')
  if (tmp(1:p-1) == ' ') return
  if (is_dble(tmp(1:p-1))) res = res + 1
  tmp = adjustlt(tmp(p+1:len_trim(tmp)))
end do
end function

!-----------------------------------------------------------------------
! lcase
!-----------------------------------------------------------------------
function lcase(str) result(res)
!! Converts a string to lower case.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: lcase`  
!! `implicit none`  
!! `print*, lcase('SoMe UpPeR')`  
!! `end program`  
!! @note This function uses alphabetic order of ASCII characters.  
character(*), intent(in) :: str !! Variable to convert to lower case.
character(len(str))      :: res !! Content of str converted to lower case.
integer :: diff, i

res = str
diff = ichar('A') - ichar('a')
do i = 1, len_trim(str)
  if (str(i:i) < 'A' .or. str(i:i) > 'Z') cycle
  res(i:i) = char(ichar(str(i:i)) - diff)
end do
end function

!-----------------------------------------------------------------------
! adjustlt
!-----------------------------------------------------------------------
function adjustlt(string) result(res)
!! Complements intrinsic procedure `adjustl` removing TAB characters.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: adjustlt`  
!! `implicit none`  
!! `print*, adjustlt('  	word')`  
!! `end program`  
!! @note This function takes the TAB character as the resutl of char(9).  
character(*), intent(in) :: string !! Character variable.
character(len(string))   :: res    !! Content of string adjusted to the left, after removing blank spaces and tab characters.

res = string
if (len_trim(string) <= 0) return
do while ((res(1:1) == char(9) .or. res(1:1) == ' ') .and. len_trim(res)>0)
  res = res(2:len_trim(res))
enddo
end function

!--------------------------------------------------------------------
! word
!--------------------------------------------------------------------
function word(str, pos) result(res)
!! Gets the `i`-th blank-separated word of a string.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: word, trim`  
!! `implicit none`  
!! `print*, trim(word('first second third, 2))`  
!! `end program`  
!!
!! @note Please use `trim(word(argument))`.  
character(*),      intent(in) :: str !! String to be analyzed.
integer, optional, intent(in) :: pos !! Position of the blank-separated word.
character(len=len(str))       :: res !! Blank-separated word.
character(len(str)) :: tmp, part
integer :: ios, n, i, lquote

n = 1
if (present(pos)) n = pos
if (n > word_count(str)) call error('(module_convers/word) requested position '//trim(string(n))//' is beyong string length: '//&
trim(string(word_count(str))))
!search n-th word
tmp = adjustlt(str)
do i = 1, n
  if (tmp(1:1) == '"' .or. tmp(1:1) == "'") then
    lquote = 2
  else
    lquote = 0
  end if
  read(tmp, *, iostat=ios) res
  if (ios /= 0) call error('(module_convers/word) unable to read from string '//trim(tmp)//', #'//trim(string(ios)))
  res = adjustlt(res)
  tmp = adjustlt(tmp)
  tmp = adjustlt(tmp(len_trim(res)+lquote+1:))
  !detect slash
  do while (tmp(1:1) == '/')
    tmp = tmp(2:)
    read(tmp, *, iostat=ios) part
    if (ios /= 0) call error('(module_convers/word) unable to read from string '//trim(tmp)//', #'//trim(string(ios)))
    part = adjustlt(part)
    tmp  = adjustlt(tmp)
    tmp  = adjustlt(tmp(len_trim(part)+1:))
    res  = trim(res)//'/'//trim(part)
  end do
  !detect comma(s)
  do while (tmp(1:1) == ',')
    tmp = adjustlt(tmp(2:))
  end do
end do
end function

!-----------------------------------------------------------------------
! trim_prv (trim)
!-----------------------------------------------------------------------
function trim_prv(str, sep, back) result(res)
!! Returns the beginning of a string before the first/last appearance of separator.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: trim`  
!! `implicit none`  
!! `print*, trim('testing.file.txt','.', back=.true.)`  
!! `end program`  
character(*),      intent(in) :: str  !! String to be analyzed.
character(*),      intent(in) :: sep  !! Separator.
logical, optional, intent(in) :: back !! If .true., search the last appearance of separator.
character(len(str)) :: res
integer :: pos

pos = index(str, sep, back = back)
if (pos == 0) then
  res = str
else
  res = str(1:pos-1)
endif
end function

!-----------------------------------------------------------------------
! replace
!-----------------------------------------------------------------------
subroutine replace(formula, cad1, cad2)
!! Replace in a string one substring by a second one (subroutine version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: reaplace`  
!! `implicit none`  
!! `character(10) :: str = 'house'`  
!! `call replace(str, 'h', 'm')`  
!! `print*, str`  
!! `end program`  
!!
!! @warning The two substrings can have different lengths. If the second substring is longer than the first one,
!! an "overflow" can occurs in the substitution.  
character(*), intent(inout) :: formula  !! String where the replacement will take place.
character(*), intent(in)    :: cad1     !! Substring to replace to.
character(*), intent(in)    :: cad2     !! Substring to replace with.
integer :: pos

pos = index(trim(formula), cad1)
do while (pos>0)
  formula = formula(1:pos-1)//cad2//formula(pos+len(cad1):len_trim(formula))
  pos = index(trim(formula), cad1)
enddo
end subroutine

!-----------------------------------------------------------------------
! freplace
!-----------------------------------------------------------------------
function freplace(formula, cad1, cad2) result(res)
!! Replace in a string one substring by a second one (function version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: freplace`  
!! `implicit none`  
!! `print*, freplace('house', 'h', 'm')`  
!! `end program`  
!!
!! @warning The two substrings can have different lengths. If the second substring is longer than the first one,
!! an "overflow" can occurs in the substitution.  
character(*), intent(in) :: formula !! String where the replacement will take place.
character(*), intent(in) :: cad1    !! Substring to replace to.
character(*), intent(in) :: cad2    !! Substring to replace with.
character(len(formula))  :: res     !! Result of the substitution.

res = formula
call replace(res, cad1,cad2)
end function

end module
