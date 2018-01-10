program test_alloc_log_r1
use module_alloc
implicit none
logical, allocatable :: v(:)
integer, allocatable :: x(:)
integer :: p

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
print*, 'call extend(v, 1, fit=.false.)'
         call extend(v, 1, fit=.false.)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, 'New extension: ', size(v,1)

print*, ' '
print*, 'Setting a scalar in a preallocated 1x3 array:'
print*, 'call alloc(v, 3)'
         call alloc(v, 3)
print*, 'call set(v, .true., 1)'
         call set(v, .true., 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Setting a vector in a preallocated 1x3 array:'
print*, 'call set(v, [.false., .true.], [1, 3])'
         call set(v, [.false., .true.], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, ' '
print*, 'Adding a scalar in a preallocated 1x3 array:'
print*, 'call add(v, .true., 1)'
         call add(v, .true., 1)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Adding a vector in a preallocated 1x3 array:'
print*, 'call add(v, [.true., .false.], [1, 3])'
         call add(v, [.true., .false.], [1, 3])
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Inserting a scalar in a preallocated 1x3 array:'
print*, 'call insert(v, .false., 5)'
         call insert(v, .false., 5)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Reducing a preallocated 1x7 array to 1x2:'
print*, 'call reduce(v, 2)'
         call reduce(v, 2)
if (.not. allocated(v)) stop 'Error, v not allocated'
print*, v

print*, 'Find the first position of a scalar value in the array:'
print*, 'call set(v,[.true.,.false.,.true.], [1,2,3]); p = find_first(v, .false.)'
         call set(v,[.true.,.false.,.true.], [1,2,3]); p = find_first(v, .false.)
print*, p

print*, 'Find all the positions of a scalar value in the array:'
print*, 'call set(v,[.true.,.false.,.true.], [1,2,3]); x = find(v)'
call dealloc(x)
         call set(v,[.true.,.false.,.true.], [1,2,3]); x = find(v) ! allocate(x, source = find(v))
if (.not. allocated(x)) stop 'Error, x not allocated'
print*, x

print*, 'Find all the positions of a scalar value in the array:'
print*, 'call set(v,[.true.,.false.,.true.], [1,2,3]); call sfind(v, .true., x)'
         call set(v,[.true.,.false.,.true.], [1,2,3]); call sfind(v, .true., x)
if (.not. allocated(x)) stop 'Error, x not allocated'
print*, x

end program   
