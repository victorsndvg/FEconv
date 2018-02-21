program test_os
use basicmod, only: set_os, get_os, slash
implicit none
character(100) :: cad

cad = get_os()
print*, cad
print*, slash()
call set_os('windows')
print*, slash()
call set_os('mac')
end program
