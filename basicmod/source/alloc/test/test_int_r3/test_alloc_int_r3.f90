program test_alloc_int_r3
use module_alloc
implicit none  
integer, allocatable :: v(:,:,:)
integer :: val(8)
integer, allocatable :: res(:,:)
call dealloc(v)
call alloc(v, 20, 30, 10)
print*, v
call extend(v, 2, 40, 20, fit=[.false.,.true., .true.])
print*, v
call reduce(v, 20, 20, 20)
print*, v
call set(v, 25, 10, 20, 5, fit=[.true.,.false.,.true.])
print*, v
call set(1, v, reshape([1,2,3,4],[2,2]), 10, fit=[.true.,.false.,.true.])
print*, v
call set(v, reshape([1,2,3,4,5,6,7,8],[2,2,2]), [10,20], [3,4], [1,3], fit=[.true.,.false.,.true.])
print*, v
call add(v, 25, 10, 20, 7, fit=[.true.,.false.,.true.])
print*, v
call add(1, v, reshape([1,2,3,4],[2,2]), 10, fit=[.true.,.false.,.true.])
print*, v
call add(v, reshape([1,2,3,4,5,6,7,8],[2,2,2]), [10,20], [3,4], [1,3], fit=[.true.,.false.,.true.])
print*, v
call insert(1, v, reshape([1,2,3,4],[2,2]), 10, used=3, fit=[.true.,.false.,.true.])
print*, v
call sfind(v, 2, res)
print*, res
call sfind(v, val, res)
print*, res
end program
