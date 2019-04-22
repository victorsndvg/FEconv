module module_alloc_bmod
!! License: GNU GPLv2
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Eases the allocation of variables.  
!!
!! @note This module does not have its own code; it only uses other modules for allocation of specific data types.  

!integer arrays
use module_alloc_int_r1_bmod
use module_alloc_int_r2_bmod
use module_alloc_int_r3_bmod

!real arrays
!use module_alloc_real_r1
!use module_alloc_real_r2
!use module_alloc_real_r2_alloc

!real64 arrays
use module_alloc_real64_r1_bmod
use module_alloc_real64_r2_bmod
use module_alloc_real64_r3_bmod

!complex64 arrays
use module_alloc_complex64_r1_bmod
use module_alloc_complex64_r2_bmod

!character arrays
use module_alloc_char_r1_bmod
use module_alloc_char_r2_bmod
!use module_alloc_char_r2_alloc

!logical arrays
use module_alloc_log_r1_bmod
!use module_alloc_log_r2
!use module_alloc_log_r2_alloc
implicit none

end module
