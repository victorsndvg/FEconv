program test_args
use module_args, only: args_count, get_arg, is_arg, get_post_arg
implicit none

! Execute in command line: test_args a -d c
print*, 'Number of arguments given in the command line: ', args_count()
print*, 'Argument n.1: ', trim(get_arg(1))
print*, 'Argument ''-d'' is given? ', is_arg('-d')
if (is_arg('-d')) then
  print*, 'Argument following ''-d'': ', trim(get_post_arg('-d'))
end if

end program   
