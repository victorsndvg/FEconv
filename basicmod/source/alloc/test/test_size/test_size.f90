program prueba
use module_alloc, only:csize
implicit none
integer :: ii1(8)
integer, allocatable :: v1(:), v2(:,:)
character(80), allocatable :: c1(:)
print*, csize(int1=ii1, d=1)
print*, csize(int1=v1,  d=1)
allocate(v2(10,10))
print*, csize(int2=v2,  d=1)
print*, csize(char1=c1, d=1)
print*, csize(d=3)
end program