module module_math_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for mathematical operations.
!
! PUBLIC PROCEDURES:
! - cross_product: calculates the cross product of 3-component vectors
! - det: find the determinant of a square matrix
! - inv2: calculates the inverse of a matrix 2x2 real64
! - inv3: calculates the inverse of a matrix 3x3 real64
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: real64
use module_report_bmod, only: error
implicit none

!Private procedures
private :: swap

contains

!-----------------------------------------------------------------------
! cross_product
!-----------------------------------------------------------------------
function cross_product(u,v) result(res)
!! Calculates the cross product of 3-component vectors.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, cross_product`  
!! `implicit none`  
!! `real(real64) :: u(3), v(3)`  
!! `u = real([1,0,0], real64)`  
!! `v = real([0,1,0], real64)`  
!! `print*, cross_product(u,v)`  
!! `end program`  
real(real64), dimension(3), intent(in) :: u   !! First vector.
real(real64), dimension(3), intent(in) :: v   !! Second vector.
real(real64), dimension(3)             :: res !! Cross product.

res = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)] 
end function

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
  p = maxloc(abs(m(k:n, k:n)), mask = abs(m(k:n, k:n))>epsilon(m))
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
! inv2
!-----------------------------------------------------------------------
subroutine inv2(A, res)
!! Calculates the inverse of a matrix 2x2 real64.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, inv2`  
!! `implicit none`  
!! `real(real64) :: A(2,2)`  
!! `real(real64) :: A(2,2), Ai(2,2)`  
!! `A = reshape(real([1,2,3,4],real64), [2,2])`  
!! `call inv2(A, Ai)`  
!! `print*, 'inv: ', Ai`  
!! `end program`  
!! @note The formula was obtained from the [Wolfram MathWorld website](http://mathworld.wolfram.com/MatrixInverse.html).  
!!
!! @warning The output matrix must be preallocated.  
real(real64), intent(in)    ::   A(2,2) !! Matrix.
real(real64), intent(inout) :: res(2,2) !! Inverse matrix.
real(real64) :: det
  
det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
if (abs(det) <= epsilon(1._real64)) call error('(module_math::inv2) unable to inverse matrix, determinant is  close to 0.')
res(1,:) = [ A(2,2), -A(1,2)]/det
res(2,:) = [-A(2,1),  A(1,1)]/det
end subroutine

!-----------------------------------------------------------------------
! inv3
!-----------------------------------------------------------------------
subroutine inv3(A, res)
!! Calculates the inverse of a matrix 3x3 real64.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, inv3`  
!! `implicit none`  
!! `real(real64) :: B(3,3), Bi(3,3)`  
!! `B = reshape(real([1,2,3,4,5,6,7,8,0],real64), [3,3])`  
!! `call inv3(B, Bi)`  
!! `print*, 'inv: ', Bi`  
!! `end program`  
!! @note The formula for the inverse was obtained from the
!! [Wolfram MathWorld website](http://mathworld.wolfram.com/MatrixInverse.html).  
!! The formula for the determinant was obtained from the
!! [Wolfram MathWorld website](https://es.mathworks.com/help/aeroblks/determinantof3x3matrix.html).  
!!
!! @warning The output matrix must be preallocated.  
real(real64), intent(in)    ::   A(3,3) !! Matrix.
real(real64), intent(inout) :: res(3,3) !! Inverse matrix.
real(real64) :: det

det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
     -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
     +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
if (abs(det) <= epsilon(1._real64)) call error('(module_math::inv3) unable to inverse matrix, determinant is  close to 0.')
res(1,:) = [A(2,2)*A(3,3)-A(3,2)*A(2,3), A(1,3)*A(3,2)-A(3,3)*A(1,2), A(1,2)*A(2,3)-A(2,2)*A(1,3)]/det
res(2,:) = [A(2,3)*A(3,1)-A(3,3)*A(2,1), A(1,1)*A(3,3)-A(3,1)*A(1,3), A(1,3)*A(2,1)-A(2,3)*A(1,1)]/det
res(3,:) = [A(2,1)*A(3,2)-A(3,1)*A(2,2), A(1,2)*A(3,1)-A(3,2)*A(1,1), A(1,1)*A(2,2)-A(2,1)*A(1,2)]/det
end subroutine
  
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
