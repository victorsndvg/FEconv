module module_math
!-----------------------------------------------------------------------
! Module for mathematical operations
! Last update: 10/04/2012
! Programmer: fran(dot)pena(at)usc(dot)es
!
! PUBLIC PROCEDURES:
! - det: find the determinant of a square matrix 
!-----------------------------------------------------------------------
use module_compiler_dependant
implicit none

!Private procedures
private :: swap

contains

!-----------------------------------------------------------------------
! det: find the determinant of a square matrix 
! The equivalent upper triangular matrix is calculated; then, det is
! the product of the diagonal
!-----------------------------------------------------------------------
function det(muser) result(res)

  real(real64), dimension(:,:), intent(in) :: muser
  real(real64), dimension(size(muser,1),size(muser,2)) :: m
  real(real64) :: res
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
      
!--------------------------------------------------------------------
! swap: swap the values of two real arrays
!--------------------------------------------------------------------
subroutine swap(u, v)

  real(real64), dimension(:), intent(inout) :: u, v
  real(real64), dimension(size(u)) :: tmp
  
  tmp = u; u = v; v = tmp
  
end subroutine

end module
