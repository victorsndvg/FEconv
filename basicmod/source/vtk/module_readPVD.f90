module module_readPVD_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es and Andres Prieto, andres.prieto(at)usc.es
!! Date: 31/12/2017
!!
!! Module to read PVD files.  
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: real64
use module_report_bmod, only: error
use module_convers_bmod, only: string
use module_files_bmod, only: get_unit
implicit none

contains

!-----------------------------------------------------------------------
! readPVD
!-----------------------------------------------------------------------
subroutine readPVD(fich, time, vtu)
!! Read a PVD file.  
!! __Example:__  
!! `program test_readPVD`  
!! `use basicmod, only: real64, maxpath, readPVD`  
!! `implicit none`  
!! `real(real64), allocatable :: time(:)`  
!! `character(maxpath), allocatable :: vtu(:)`  
!! `call readPVD('temp_evolution.pvd', time, vtu)`  
!! `end program`  
character(*),               intent(in)  :: fich    !! PVD filename.
real(real64),  allocatable, intent(out) :: time(:) !! Time array.
character(*),  allocatable, intent(out) :: vtu(:)  !! VTU filename array.
character(255) :: line
integer :: un, ios, p, q, n, i

un = get_unit()
open (unit=un, file=fich, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('(module_readPVD/readPVD) open error #'//trim(string(ios)))
n = 0
do
  read (unit=un, fmt='(a)', iostat=ios) line
  if (ios < 0) exit !found EOF
  if (index(line, 'timestep="') > 0) n = n + 1
end do
allocate(time(n), vtu(n))
rewind(un)
i = 0
do
  if (i == n) exit
  read (unit=un, fmt='(a)', iostat=ios) line
  if (ios < 0) call error('readPVD/read, #'//trim(string(ios)))
  !time
  p = index(line, 'timestep="')
  if (p > 0) then
    i = i + 1
    q = index(line(p+10:),'"')
    read(line(p+10:q+p+8),*) time(i)
    !file
    p = index(line, 'file="')
    if (p == 0) stop 'ERROR: mark file= not found'
    q = index(line(p+6:),'"')
    vtu(i) = line(p+6:q+p+4)
  end if
end do
end subroutine

end module
