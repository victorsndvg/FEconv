program test_writevtu
use module_writePVD, only: real64, writePVD
implicit none
real(real64) :: time(5) = [0._real64, 1._real64, 1.5_real64, 10._real64, 11.5_real64]
call writePVD('temp_evolution.pvd', 'temperature', time)
end program