program test
use module_compiler_dependant
implicit none
integer :: exits, cmds
call execute_command_line('program.exe', exitstat=exits, cmdstat=cmds)
if (cmds /= 0) write(error_unit,'(a)') cmds
print*, exits
end program
