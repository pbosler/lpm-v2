module STDIntVectorModule
!> @file STDIntVector.f90
!> Provides a class and methods that mimic a C++ `std::vector<int>` object.  
!> 
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup STDIntVector STDIntVector module
!> @brief Mimics the dynamically allocated C++ `std::vector<int>` object.  
!> Provides array-style access to its elements, and resizes if necessary.  
!> @{
use NumberKindsModule

implicit none
private
public STDIntVector, initialize

type STDIntVector
	integer(kint), pointer :: integers(:) => null()
	integer(kint) :: N = 0
	
	contains
		procedure :: defaultInit
		procedure :: initFromArray
		procedure :: initFromCopy
		final :: finalizeVector
		procedure, public :: empty
		procedure, public :: int
		procedure, public :: pushBack
		procedure, public :: pushBackUnique
		procedure, public :: insert
		procedure, public :: replace
		procedure, public :: print
end type

integer(kint), parameter :: DEFAULT_VEC_SIZE = 20

interface initialize
	module procedure defaultInit
	module procedure initFromArray
	module procedure initFromCopy
end interface

contains

subroutine defaultInit(self, minSize)
	class(STDIntVector), intent(out) :: self
	integer(kint), intent(in), optional :: minSize
	if ( present(minSize) .AND. minSize > DEFAULT_VEC_SIZE ) then
		allocate(self%integers(minSize))
	else
		allocate(self%integers(DEFAULT_VEC_SIZE))
	endif
	self%integers = 0
	self%N = 0
end subroutine

subroutine initFromArray( self, intArray ) 
	class(STDIntVector) :: self
	integer(kint), intent(in) :: intArray(:)
	
	if ( size(intArray) < DEFAULT_VEC_SIZE ) then
		allocate(self%integers(DEFAULT_VEC_SIZE))
	else
		allocate(self%integers(size(intArray) + DEFAULT_VEC_SIZE))
	endif
	
	self%integers(1:size(intArray)) = intArray
	self%N = size(intArray)
	self%integers(size(intArray)+1:size(self%integers)) = 0	
end subroutine

subroutine initFromCopy( self, other ) 
	class(STDIntVector), intent(out) :: self
	class(STDIntVector), intent(in) :: other
	allocate(self%integers(size(other%integers)))
	self%N = other%N
	self%integers = other%integers
end subroutine

subroutine assignVectorFromCopy( self, other )
	type(STDIntVector), intent(out) :: self
	type(STDIntVector), intent(in) :: other
	call initFromCopy(self, other)
end subroutine

subroutine finalizeVector( vec )
	type(STDIntVector), intent(inout) :: vec
	if ( associated(vec%integers)) deallocate(vec%integers)
end subroutine

function empty( self ) result(tf)
	logical(klog) :: tf
	class(STDIntVector), intent(in) :: self
	tf = ( self%N == 0 )
end function

function int(self, index) result(j)
	integer(kint) :: j
	class(STDIntVector), intent(in) :: self
	integer(kint), intent(in) :: index
	j = 0
	if ( index < 1 .OR. index > self%N) then
		print *,"STDIntVector%int ERROR : attempted to access index out of bounds."
	else
		j = self%integers(index)
	endif
end function

subroutine pushBack(self, int)
	class(STDIntVector), intent(inout) :: self
	integer(kint), intent(in) :: int
	
	if ( self%N+1 > size(self%integers) ) &
		call doubleArraySpace(self)
	self%integers(self%N+1) = int
	self%N = self%N + 1
end subroutine

subroutine pushBackUnique( self, int )
	class(STDIntVector), intent(inout) :: self
	integer(kint), intent(in) :: int
	integer(kint) :: i
	logical(klog) :: isNew
	
	isNew = .TRUE.
	do i = 1, self%N
		if ( int == self%integers(i) ) isNew = .FALSE.
	enddo
	
	if ( isNew ) call pushBack(self, int)
end subroutine

subroutine insert( self, pos, val )
	class(STDIntVector), intent(inout) :: self
	integer(kint), intent(in) :: pos
	integer(kint), intent(in) :: val
	!
	integer(kint), allocatable :: tempArray(:)
	
	if ( pos < 0 .OR. pos > self%N+1 ) then
		print *, "STDIntVector%insert ERROR : attempted to access index out of bounds."
		return
	endif
	
	if ( self%N + 1 > size(self%integers) ) &
		call doubleArraySpace(self)
		
	if ( pos == self%N + 1 ) then
		call pushBack(self, val)
	else	
		allocate(tempArray(self%N))
		tempArray = self%integers(1:self%N)
	
		self%integers(1:pos-1) = tempArray(1:pos-1)
		self%integers(pos) = val
		self%integers(pos+1:self%N+1) = tempArray(pos:self%N)
	
		self%N = self%N + 1	
		deallocate(tempArray)
	endif
end subroutine

subroutine replace( self, pos, val )
	class(STDIntVector), intent(inout) :: self
	integer(kint), intent(in) :: pos
	integer(kint), intent(in) :: val
	
	if ( pos < 0 .OR. pos > self%N ) then
		print *, "STDIntVector%replace ERROR : attempted to access index out of bounds."
		return
	endif
	self%integers(pos) = val	
end subroutine

subroutine doubleArraySpace( self )
	class(STDIntVector), intent(inout) :: self
	!
	integer(kint), allocatable :: tempArray(:)
	integer(kint) :: m
	
	m = size(self%integers)
	allocate(tempArray(self%N))
	tempArray = self%integers(1:self%N)
	
	deallocate(self%integers)
	allocate(self%integers( 2 * m ))
	self%integers(1:self%N) = tempArray
	self%integers(self%N+1:2*m) = 0	
	deallocate(tempArray)
end subroutine

subroutine print(self, name )
	class(STDIntVector), intent(in) :: self
	character(len=*), intent(in), optional :: name
	!
	integer(kint) :: i
	
	if ( present(name) ) then
		write(6,'(A,I6)') trim(name)//"%N = ", self%N
	endif
	do i = 1, self%N
		write(6,'(I8)',advance='NO') self%integers(i)
	enddo
	write(6,'(A)',advance='YES') ' '
end subroutine

!> @}
end module