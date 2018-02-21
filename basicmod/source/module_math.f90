module module_math_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for mathematical operations.
!
! PUBLIC PROCEDURES:
! - det: find the determinant of a square matrix
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: real64
implicit none

!Private procedures
private :: swap

contains

!-----------------------------------------------------------------------
! det
!-----------------------------------------------------------------------
function det(muser) result(res)
!! Calculates the determinant of a square matrix.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, det`  
!! `implicit none`  
!! `real(real64) :: A(2,2)`  
!! `A = reshape([1._real64, 2._real64, 3._real64, 4._real64], [2,2])`  
!! `print*, det(A)`  
!! `end program`  
!!
!! @note The equivalent upper triangular matrix is calculated; the determinant is just the product of the diagonal.
real(real64), dimension(:,:), intent(in) :: muser !! Square matrix.
real(real64)                             :: res   !! Determinant.
real(real64), dimension(size(muser,1),size(muser,2)) :: m
integer :: sgn, n, k, j
integer, dimension(2) :: p

res = 0; sgn = 1; m = muser; n = size(m, 1)
do k = 1, n-1
  p = maxloc(m(k:n, k:n), mask = abs(m(k:n, k:n))>epsilon(m))
  if (all(p == 0)) return !no valid pivot
  p = p + k - 1 !position respect to the total matrix
  if (p(1) /= k) then !permute rows
    call swap(m(k,:), m(p(1),:)); sgn = -sgn
  end if
  if (p(2) /= k) then !permute cols
    call swap(m(:,k), m(:,p(2))); sgn = -sgn
  end if
  do j = k+1, n !use the pivot
    m(j,:) = m(j,:) - m(j,k) * m(k,:) / m(k,k)
  end do
end do
res = sgn * product((/(m(k,k), k=1,n)/)) !calculation of det
end function

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!--------------------------------------------------------------------
! swap
!--------------------------------------------------------------------
subroutine swap(u, v)
!! Swaps the values of two real arrays.
real(real64), dimension(:), intent(inout) :: u, v
real(real64), dimension(size(u)) :: tmp

tmp = u; u = v; v = tmp
end subroutine

end module
