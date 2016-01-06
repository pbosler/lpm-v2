module STDIntVectorModule
!> @file STDIntVector.f90
!> Provides a class and methods that mimic a C++ `std::vector<int>` object.  
!> 
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup STDIntVector STDIntVector
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

integer(kint), parameter :: DEFAULT_VEC_SIZE = 20 !< default size of array

!< @brief Interface to initialize a STDIntVector
interface initialize
	module procedure defaultInit
	module procedure initFromArray
	module procedure initFromCopy
end interface

contains

!< @brief Initializes an empty STDIntVector.  Size and vector entries are zeroed.
!> @param self Target to be initialized
!> @param minSize minimum size of vector
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

!< @brief Initializes an STDIntVector with the contents and size of an integer array.
!> @param self Target to be initialized
!> @param intArray integer array with source data
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

!< @brief Initializes an STDIntVector by copying another STDIntVector
!> @param self Target to be initialized
!> @param other STDIntVector with source data
subroutine initFromCopy( self, other ) 
	class(STDIntVector), intent(out) :: self
	class(STDIntVector), intent(in) :: other
	allocate(self%integers(size(other%integers)))
	self%N = other%N
	self%integers = other%integers
end subroutine

!< @brief Assignment operator function 
!> @todo Implement assignment operator
!> 
!> @param self Target to be overwritten with content of other
!> @param other Source STDIntvector
subroutine assignVectorFromCopy( self, other )
	type(STDIntVector), intent(out) :: self
	type(STDIntVector), intent(in) :: other
	call initFromCopy(self, other)
end subroutine

!> @brief Deletes and frees memory associated with a STDIntVector object.
!> @param vec Target STDIntVector
subroutine finalizeVector( vec )
	type(STDIntVector), intent(inout) :: vec
	if ( associated(vec%integers)) deallocate(vec%integers)
end subroutine

!> @brief Returns .TRUE. if a STDIntVector has no content.
!> @param self Target STDIntVector
!> @return .TRUE. if self is empty
function empty( self ) result(tf)
	logical(klog) :: tf
	class(STDIntVector), intent(in) :: self
	tf = ( self%N == 0 )
end function

!> @brief Returns the integer in the STDIntVector at the desired position
!> @param self Source STDIntVector
!> @param index Location of desired integer in STDIntVector's array
!> @return The integer in the array at the index location
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

!> @brief Adds a new integer to the end of a vector, resizes if necessary.
!> @param self STDIntVector to increase in size
!> @param int Integer to add to self
subroutine pushBack(self, int)
	class(STDIntVector), intent(inout) :: self
	integer(kint), intent(in) :: int
	
	if ( self%N+1 > size(self%integers) ) &
		call doubleArraySpace(self)
	self%integers(self%N+1) = int
	self%N = self%N + 1
end subroutine

!> @brief Adds a new integer to the end of a vector, but only if that integer is not already listed somewhere in the vector.
!> @param self Target STDIntVector
!> @param int  Integer to add (if it is unique)
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

!> @brief Adds a new integer to a STDIntVector at a desired position; increases the size of the vector.
!> @param self Target STDIntVector
!> @param pos Position of new integer in array
!> @param val Integer to add at new position
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

!> @brief Replaces the value of an integer at a particular position in the array with a new value.  
!> Does not change the size of the STDIntVector
!> @param self Target STDIntVector
!> @param pos Position of integer to replace
!> @param val New integer to put into array
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

!> @brief Doubles the memory assigned to a STDIntVector.
!> @param self Target STDIntVector
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

!< @brief Prints the contents of a vector to the console.
!> @param self Target STDIntVector
!> @param name Name of vector
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