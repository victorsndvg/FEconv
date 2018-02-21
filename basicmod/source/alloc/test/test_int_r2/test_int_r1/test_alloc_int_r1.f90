program test_alloc_int_r1
use module_alloc
implicit none
integer, allocatable :: v(:)

print*, 'Extension of a non-allocated array without fitting:'
print*, 'call extend(v, 1, fit=.false.)'
         call extend(v, 1, fit=.false.)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extension: ', size(v,1)
deallocate(v)

print*, ' '
print*, 'Extension of a non-allocated array with fitting:'
print*, 'call extend(v, 1)'
         call extend(v, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extension: ', size(v,1)

print*, ' '
print*, 'Extension of an allocated array with fitting:'
print*, 'call extend(v, 2)'
         call extend(v, 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extension: ', size(v,1)

print*, ' '
print*, 'Extension of an allocated array without surpassing current dimensions:'
print*, 'call extend(v, 1, .false.)'
         call extend(v, 1, .false.)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extension: ', size(v,1)

print*, ' '
print*, 'Setting a scalar in a preallocated 1x3 array:'
print*, 'call alloc(v, 3)'
         call alloc(v, 3)
print*, 'call set(v, 23, 1)'
         call set(v, 23, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Setting a vector in a preallocated 1x3 array:'
print*, 'call set(v, [4, 5], [1, 3])'
         call set(v, [4, 5], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, ' '
print*, 'Adding a scalar in a preallocated 1x3 array:'
print*, 'call add(v, 23, 1)'
         call add(v, 23, 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Adding a vector in a preallocated 1x3 array:'
print*, 'call add(v, [4, 5], [1, 3])'
         call add(v, [4, 5], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Inserting a scalar in a preallocated 1x3 array:'
print*, 'call insert(v, 41, 5)'
         call insert(v, 41, 5)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Reducing a preallocated 1x7 array to 1x2:'
print*, 'call reduce(v, 2)'
         call reduce(v, 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

end program   
