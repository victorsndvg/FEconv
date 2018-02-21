program test_xread
use module_report, only: report_option
use module_xread
implicit none
integer :: res, a
integer, allocatable :: i(:)
real(real64) :: x
real(real64), allocatable :: x2(:)
character(maxpath) :: c
character(maxpath), allocatable :: c2(:), val(:)

call report_option('output', 'std')
print*, 'Reading data.xml...'
open(10, file='data.xml')
res = xread(10, 'cosa',     a,  rew=.true.)
res = xread(10, 'dato',     i,  rew=.true.)
res = xread(10, 'constant', x,  rew=.true.)
res = xread(10, 'c',        x2, rew=.true.)
res = xread(10, 'mesh',     c,  rew=.true.)
res = xread(10, 'm',        c2, rew=.true.)
close(10)

print*, 'Reading data.xml...'
open(10, file='data.xml')
call  xlist(10, '.data', val, rew=.true.)
res = xread(10, '.data.'//val(1), a,  rew=.true.)
res = xread(10, '.data.'//val(2), i,  rew=.true.)
res = xread(10, '.data.'//val(3), x,  rew=.true.)
res = xread(10, '.data.'//val(4), x2, rew=.true.)
res = xread(10, '.data.'//val(5), c,  rew=.true.)
res = xread(10, '.data.'//val(6), c2, rew=.true.)
close(10)

end program
