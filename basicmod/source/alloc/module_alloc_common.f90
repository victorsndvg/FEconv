module module_alloc_common_bmod
!! License: GNU GPLv2
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module that contains common entities used by all allocation modules.  
use module_compiler_dependant_bmod, only: real64
implicit none

!Constants
integer, parameter :: DEFAULT_ALLOC  = 1000 !! Initial size for allocation.

contains

!-----------------------------------------------------------------------
! csize
!-----------------------------------------------------------------------
function csize(int1, int2, int3, real1, real2, real3, dbl1, dbl2, dbl3, cmplx1, cmplx2, cmplx3, char1, char2, char3, log1, log2, &
log3, d) result(res)
!! Returns zero for unallocated arrays and their size otherwise.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: csize`  
!! `implicit none`  
!! `integer :: x(5)`  
!! `real, allocatable :: y(:), z(:,:)`  
!! `allocate(y(10))`  
!! `print*, csize(int1=x,  d=1) ! prints  5`  
!! `print*, csize(real1=y, d=1) ! prints 10`  
!! `print*, csize(real2=z, d=1) ! prints  0`  
!! `end program`  
!!
!! @note Fortran standard states that an unallocated allocatable does not have size and thus it is illegal to check it. 
!! Procedure `csize` returns zero in such a case. This is achieved using the Fortan 2008 enhancement that states that a 
!! null pointer or an unallocated allocatable can be used to denote an absent non-allocatable non-pointer optional argument.  
!!
!! @warning Since the array arguments are optional, you must use explicit keywords to avoid type mismatch.  
integer,         intent(in), optional :: int1(:)       !! Integer   rank 1 array.
integer,         intent(in), optional :: int2(:,:)     !! Integer   rank 2 array.
integer,         intent(in), optional :: int3(:,:,:)   !! Integer   rank 3 array.
real,            intent(in), optional :: real1(:)      !! Real      rank 1 array.
real,            intent(in), optional :: real2(:,:)    !! Real      rank 2 array.
real,            intent(in), optional :: real3(:,:,:)  !! Real      rank 3 array.
real(real64),    intent(in), optional :: dbl1(:)       !! Real64    rank 1 array.
real(real64),    intent(in), optional :: dbl2(:,:)     !! Real64    rank 2 array.
real(real64),    intent(in), optional :: dbl3(:,:,:)   !! Real64    rank 3 array.
complex(real64), intent(in), optional :: cmplx1(:)     !! Complex64 rank 1 array.
complex(real64), intent(in), optional :: cmplx2(:,:)   !! Complex64 rank 2 array.
complex(real64), intent(in), optional :: cmplx3(:,:,:) !! Complex64 rank 3 array.
character(*),    intent(in), optional :: char1(:)      !! Character rank 1 array.
character(*),    intent(in), optional :: char2(:,:)    !! Character rank 2 array.
character(*),    intent(in), optional :: char3(:,:,:)  !! Character rank 3 array.
logical,         intent(in), optional :: log1(:)       !! Logical   rank 1 array.
logical,         intent(in), optional :: log2(:,:)     !! Logical   rank 2 array.
logical,         intent(in), optional :: log3(:,:,:)   !! Logical   rank 3 array.
integer,         intent(in)           :: d             !! Dimension
integer                               :: res           !! 0 for an unallocated arrays; their size otherwise.

if     (present(int1))   then; res = size(int1,   d)
elseif (present(int2))   then; res = size(int2,   d)
elseif (present(int3))   then; res = size(int3,   d)
elseif (present(real1))  then; res = size(real1,  d)
elseif (present(real2))  then; res = size(real2,  d)
elseif (present(real3))  then; res = size(real3,  d)
elseif (present(dbl1))   then; res = size(dbl1,   d)
elseif (present(dbl2))   then; res = size(dbl2,   d)
elseif (present(dbl3))   then; res = size(dbl3,   d)
elseif (present(cmplx1)) then; res = size(cmplx1, d)
elseif (present(cmplx2)) then; res = size(cmplx2, d)
elseif (present(cmplx3)) then; res = size(cmplx3, d)
elseif (present(char1))  then; res = size(char1,  d)
elseif (present(char2))  then; res = size(char2,  d)
elseif (present(char3))  then; res = size(char3,  d)
elseif (present(log1))   then; res = size(log1,   d)
elseif (present(log2))   then; res = size(log2,   d)
elseif (present(log3))   then; res = size(log3,   d)
else; res = 0
end if
end function

!-----------------------------------------------------------------------
! search_multiple
!-----------------------------------------------------------------------
integer function search_multiple(a,b)
  !! Searchs the smallest value of 2 to the power of a that is bigger than b 
  !! 2^n*a > b  <=>  n > log2(b/a).
integer, intent(in) :: a, b

if (b > a) then
  search_multiple = int(2**real(ceiling(log(real(b)/a)/log(2.)))*a)
else
  search_multiple = a
end if
end function

end module
