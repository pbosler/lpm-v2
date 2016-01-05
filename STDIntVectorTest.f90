program stdIntVectorTest

use NumberKindsModule
use STDIntVectorModule

implicit none

type(STDIntVector) :: intVec1, intVec2, intVec3, intVec4

integer(kint) :: i, j
integer(kint), dimension(:) :: intArray(10), tooBigIntArray(30)


print *, "STDIntVectorTest: 1. test all constructors."
print *, "STDIntVectorTest: 2. test copy."
print *, "STDIntVectorTest: 3. test pushback, insert, replace."
print *, "STDIntVectorTest: 4. test resize."

!
! constructor tests
!
call initialize(intVec1)
print *, "Default contstructor : "
if ( .NOT. associated(intVec1%integers) ) then
	print *, "ERROR : default constructor."
endif
print *, "size intVec1%integers = ", size(intVec1%integers)
print *, "intVec1%N = ", intVec1%N
if ( intVec1%empty() ) then
	print *, "intVec1.empty() = .TRUE."
endif

print *, " "

call initialize(intVec2, 22)
if ( .NOT. associated(intVec2%integers) ) then
	print *, "ERROR : default constructor with optional arg."
endif
print *, "size intVec2%integers = ", size(intVec2%integers)
print *, "intVec2%N = ", intVec2%N

print *, " "

do i = 1, 10
	intArray(i) = 10 * i
	call intVec1%pushBack(i)
enddo

do i = 1, 30
	tooBigIntArray(i) = - i
enddo


print *, "contents of intVec1:"
do i = 1, 10
	write(6,'(I6)',advance='NO') intVec1%int(i)
enddo
print *, " "
print *, " "

print *, "inserting a number to intVec1: "
call intVec1%insert(5, 99)
print *, "contents of intVec1:"
do i = 1, intVec1%N
	write(6,'(I6)',advance='NO') intVec1%int(i)
enddo
print *, " " 
print *, " "

call initialize(intVec3, intArray)
if ( .NOT. associated(intVec3%integers) ) then
	print *, "ERROR : array constructor."
endif
print *, "size intVec3%integers = ", size(intVec3%integers)
print *, "intVec3%N = ", intVec3%N
call intVec3%replace(5,-99)
do i = 1, 10
	write(6,'(I6)',advance='NO') intVec3%int(i)
enddo
print *, " "
print *, " "

print *, "testing pushback and resize."
call initialize(intVec1)
do i = 1, 30
	call intVec1%pushback(tooBigIntArray(i))
enddo
print *, "contents of intVec1:"
do i = 1, intVec1%N
	write(6,'(I4)',advance='NO') intVec1%int(i)
enddo
print *, " "
print *, " "

print *, "testing copy constructor."
call initialize(intVec4, intVec1)
do i = 1, intVec4%N
	write(6,'(I4)',advance='NO') intVec4%int(i)
enddo
print *, " "
print *, " "




end program