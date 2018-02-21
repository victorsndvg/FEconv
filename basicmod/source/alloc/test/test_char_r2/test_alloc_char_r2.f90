program test_alloc_char_r2
use module_alloc
implicit none
character(20), allocatable :: v(:,:) 
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
         call set(v, '23', 1, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Setting a row in a preallocated 2x3 array:'
print*, 'call set(1, v, [2, 3], 2)'
         call set(1, v, ['2', '3'], 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Setting a col in a preallocated 2x3 array:'
print*, 'call set(2, v, [21, 31], 3)'
         call set(2, v, ['21', '31'], 3)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Setting a matrix in a preallocated 2x3 array:'
print*, 'call set(v, reshape([4, 5, 6, 7],[2,2]), [1, 2], [1, 3])'
         call set(v, reshape(['4', '5', '6', '7'],[2,2]), [1, 2], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, ' '
print*, 'Inserting a row in a preallocated 2x3 array:'
print*, 'call insert(1, v, [40, 50], 2)'
         call insert(1, v, ['40', '50'], 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
do i = 1, size(v,1)
  print*, (v(i,j), j = 1, size(v,2))
end do  

print*, 'Inserting a col in a preallocated 3x3 array:'
print*, 'call insert(2, v, [41, 51], 5)'
         call insert(2, v, ['41', '51'], 5)
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

end program   
