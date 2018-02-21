program test_alloc_real64_r2
use module_alloc
implicit none
real(real64), allocatable :: v(:,:) 
integer, allocatable :: res(:,:)
integer :: i, j, n

print*, 'Extension of a non-allocated array without fitting:'
print*, 'call extend(v, 1, 1, fit=[.false., .false.])'
         call extend(v, 1, 1, fit=[.false., .false.])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extensions: ', size(v,1), size(v,2)
deallocate(v)

print*, ' '
print*, 'Extension of a non-allocated array with partial fitting:'
print*, 'call extend(v, 1, 1, [.true., .false.])'
         call extend(v, 1, 1, [.true., .false.])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extensions: ', size(v,1), size(v,2)
deallocate(v)

print*, ' '
print*, 'Extension of a non-allocated array with fitting:'
print*, 'call extend(v, 1, 1)'
         call extend(v, 1, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extensions: ', size(v,1), size(v,2)

print*, ' '
print*, 'Extension of an allocated array with fitting:'
print*, 'call extend(v, 2, 3)'
         call extend(v, 2, 3)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extensions: ', size(v,1), size(v,2)

print*, ' '
print*, 'Extension of an allocated array with partial fitting:'
print*, 'call extend(v, 10, 10, [.true., .false.])'
         call extend(v, 10, 10, [.true., .false.])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extensions: ', size(v,1), size(v,2)

print*, ' '
print*, 'Extension of an allocated array without surpassing current dimensions:'
print*, 'call extend(v, 1, 1, [.false., .false.])'
         call extend(v, 1, 1, [.false., .false.])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extensions: ', size(v,1), size(v,2)

print*, ' '
print*, 'Setting a scalar in a preallocated 2x3 array:'
print*, 'call alloc(v, 2, 3)'
         call alloc(v, 2, 3)
print*, 'call set(v, 23, 1, 1)'
         call set(v, 23._real64, 1, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Setting a row in a preallocated 2x3 array:'
print*, 'call set(1, v, [2, 3], 2)'
         call set(1, v, [2._8, 3._8], 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Setting a col in a preallocated 2x3 array:'
print*, 'call set(2, v, [21, 31], 3)'
         call set(2, v, [21._8, 31._8], 3)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Setting a matrix in a preallocated 2x3 array:'
print*, 'call set(v, reshape([4, 5, 6, 7],[2,2]), [1, 2], [1, 3])'
         call set(v, reshape([4._8, 5._8, 6._8, 7._8],[2,2]), [1, 2], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, ' '
print*, 'Adding a scalar in a preallocated 2x3 array:'
print*, 'call add(v, 23, 1, 1)'
         call add(v, 23._8, 1, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Adding a row in a preallocated 2x3 array:'
print*, 'call add(1v, [2, 3], 2)'
         call add(1, v, [2._8, 3._8], 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Adding a col in a preallocated 2x3 array:'
print*, 'call add(2, v, [21, 31], 3)'
         call add(2, v, [21._8, 31._8], 3)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Adding a matrix in a preallocated 2x3 array:'
print*, 'call add(v, reshape([4, 5, 6, 7],[2,2]), [1, 2], [1, 3])'
         call add(v, reshape([4._8, 5._8, 6._8, 7._8],[2,2]), [1, 2], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, ' '
print*, 'Inserting a row in a preallocated 2x3 array:'
print*, 'call insert(1, v, [40, 50], 2)'
         call insert(1, v, [40._8, 50._8], 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Inserting a col in a preallocated 3x3 array:'
print*, 'call insert(2, v, [41, 51], 5)'
         call insert(2, v, [41._8, 51._8], 5)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Reducing a preallocated 3x5 array to 2x2:'
print*, 'call reduce(v, 2, 2)'
         call reduce(v, 2, 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Insert a row in a row-sorted array:'
print*, 'call extend(v, 2, 3); print*, size(v)'
         call extend(v, 2, 3); print*, size(v,1), size(v,2)
print*, 'call set(1, v, [1, 2, 1], 1); call set(1, v, [3, 10, 1], 2); print*, v'
         call set(1, v, [1._8, 2._8, 1._8], 1); call set(1, v, [3._8, 10._8, 1._8], 2)
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  
print*, 'n = 2; call insert_sorted(1, v, [2,20,200], used=n); print*, n, v'
         n = 2; call insert_sorted(1, v, [2._8,20._8,200._8], used=n); print*, n
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Find the position of a row in a row-sorted array:'
print*, 'print*, find_sorted(1, v, [2,20,200], used=3)'
         print*, find_sorted(1, v, [2._8,20._8,200._8], used=3)

print*, 'Find all the occurrences of scalar value in an array:'
print*, 'call sfind(v, 2, res); print*,res'
         call sfind(v, 2._8, res)
do i = 1, size(res,1)
  print*, (res(i,j), j = 1, size(res,2))
end do  

print*, 'Find all the occurrences of an array of values in an array:'
print*, 'call sfind(v, [1,2], res); print*,res'
         call sfind(v, [1._8,2._8], res)
do i = 1, size(res,1)
  print*, (res(i,j), j = 1, size(res,2))
end do  

end program   
