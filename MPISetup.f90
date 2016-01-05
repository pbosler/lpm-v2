module MPISetupModule

use NumberKindsModule

implicit none

private
public MPISetup, New, Delete, Copy
public LoadBalance

type MPISetup
	integer(kint), allocatable :: indexStart(:)
	integer(kint), allocatable :: indexEnd(:)
	integer(kint), allocatable :: messageLength(:)
	integer(kint) :: n = 0
	contains
		final :: deletePrivate
end type

interface New
	module procedure newPrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

contains

subroutine newPrivate( self, nItems, nProcs )
	type(MPISetup), intent(out) :: self
	integer(kint), intent(in) :: nItems
	integer(kint), intent(in) :: nProcs
	
	allocate(self%indexStart(0:nProcs-1))
	allocate(self%indexEnd(0:nProcs-1))
	allocate(self%messageLength(0:nProcs-1))
	
	call LoadBalance(self, nItems, nProcs)	
end subroutine

subroutine deletePrivate(self)
	type(MPISetup), intent(inout) :: self
	if ( allocated(self%indexStart) ) then
		deallocate(self%indexStart)
	 	deallocate(self%indexEnd)
	 	deallocate(self%messageLength)
	endif
end subroutine

subroutine copyPrivate(self, other)
	type(MPISetup), intent(inout) :: self
	type(MPISetup), intent(in) :: other
	call deletePrivate(self)
	call newPrivate(self, other%n, numProcs)
end subroutine

subroutine LoadBalance( self, nItems, nProcs )
	type(MPISetup), intent(inout) :: self
	integer(kint), intent(in) :: nItems, nProcs
	!
	integer(kint) :: i, chunkSize
	
	chunkSize = nItems / nProcs
	do i = 0, nProcs - 1
		self%indexStart(i) = i * chunkSize + 1
		self%indexEnd(i) = (i+1) * chunkSize
	enddo
	self%indexEnd(nProcs-1) = nItems
	self%messageLength = self%indexEnd - self%indexStart + 1
	self%n = nItems
end subroutine

end module