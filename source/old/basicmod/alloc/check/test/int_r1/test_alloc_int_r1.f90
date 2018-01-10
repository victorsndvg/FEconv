program test_alloc_int_r1
use module_alloc
implicit none
integer, allocatable :: v(:)

print*, 'Deallocation of a non-allocated array:'
print*, 'call dealloc(v)'
         call dealloc(v)

print*, ' '
print*, 'Allocation of an allocated array:'
print*, 'allocate(v(5))'
         allocate(v(5))
print*, 'call alloc(v, 1)'
         call alloc(v, 1)

print*, ' '
print*, 'Setting a scalar in a preallocated 1x5 array:'
print*, 'call set(v, 0, 0)'
         call set(v, 0, 0)
print*, 'call set(v, 23, 7)'
         call set(v, 23, 7)

print*, ' '
print*, 'Setting a vector in a preallocated 1x5 array:'
print*, 'call set(v, [4, 5], [1, 8])'
         call set(v, [4, 5], [1, 8])
print*, 'call set(v, [4, 5], [0, 8])'
!        call set(v, [4, 5], [0, 8])
print*, '(RESULT) INFO: (module_alloc_int_r1/set_vector) Given dimensions: [0, 8], size(v): 8'
print*, 'call set(v, [4, 5], [1])'
!        call set(v, [4, 5], [1])
print*, '(RESULT) INFO: (module_alloc_int_r1/set_vector) Given index dimension: 1, size(val): 2'
print*, '(RESULT) Violación de segmento'

print*, ' '
print*, 'Adding a scalar in a preallocated 1x5 array:'
         call alloc(v,5)
print*, 'call add(v, 0, 0)'
         call add(v, 0, 0)
print*, 'call add(v, 23, 7)'
         call add(v, 23, 7)

print*, ' '
print*, 'Adding a vector in a preallocated 1x5 array:'
print*, 'call add(v, [4, 5], [1, 8])'
         call add(v, [4, 5], [1, 8])
print*, 'call add(v, [4, 5], [0, 8])'
!         call add(v, [4, 5], [0, 8])
print*, '(RESULT) INFO: (module_alloc_int_r1/add_vector) Given dimensions: [0, 8], size(v): 8'
print*, 'call add(v, [4, 5], [1])'
         call add(v, [4, 5], [1])
print*, '(RESULT) INFO: (module_alloc_int_r1/add_vector) Given index dimension: 1, size(val): 2'
print*, '(RESULT) Violación de segmento'


print*, ' '
print*, 'Inserting a scalar in a preallocated array:'
print*, 'call insert(v, 23, 10)'
         call insert(v, 23, 10)
print*, 'call insert(v, 23, 0)'
         call insert(v, 23, 0)
print*, '(RESULT) INFO: (module_alloc_int_r1/insert) Given dimension: 0, size(v): 10'

end program   
