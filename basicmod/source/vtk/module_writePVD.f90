module module_writePVD_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es and Andres Prieto, andres.prieto(at)usc.es
!! Date: 31/12/2017
!!
!! Module for writing PVD files.  
!!
!! @note `writePVD` annotates every pair (_`time`(i)_, *`prefix`_i.vtu*), where *i = 1, size(`time`,1)*.  
use module_compiler_dependant_bmod, only: real64
use module_report_bmod, only: error
use module_convers_bmod, only: string
use module_files_bmod, only: get_unit

implicit none

contains

!-----------------------------------------------------------------------
! writePVD
!-----------------------------------------------------------------------
subroutine writePVD(file, prefix, time)
!! Write a PVD file.  
!! __Example:__  
!! `program test_writevtu`  
!! `use basicmod, only: real64, writePVD`  
!! `implicit none`  
!! `real(real64) :: time(5) = [0._real64, 1._real64, 1.5_real64, 10._real64, 11.5_real64]`  
!! `call writePVD('temp_evolution.pvd', 'temperature', time)`  
!! `end program`  
!!
!! @note `writePVD` annotates every pair (_`time`(i)_, *`prefix`_i.vtu*), where *i = 1, size(`time`,1)*.  
character(*), intent(in) :: file    !! PVD filename
character(*), intent(in) :: prefix  !! VTU prefix name; expected VTU files are *`prefix`_i.vtu*.
real(real64)             :: time(:) !! Time array.
integer :: un, ios, i

un = get_unit()
open (unit=un, file=file, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('(module_writePVD/writePVD) open error #'//trim(string(ios)))
!HEAD
write(un,'(a)') '<?xml version="1.0"?>'
write(un,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
write(un,'(a)') '  <Collection>'
!TIME LOOP
do i = 1,size(time,1)
  write(un,'(a)') '    <DataSet timestep="'//trim(string(real(time(i))))//'" group="" part="0" file="'//&
  &trim(prefix)//'_'//trim(string(i))//'.vtu"/>'
enddo
!TAIL
write(un,'(a)') '  </Collection>'
write(un,'(a)') '</VTKFile>'
close(un)
end subroutine

end module
