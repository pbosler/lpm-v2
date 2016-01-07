module MPISetupModule
!------------------------------------------------------------------------------
! Lagrangian Particle Method
!------------------------------------------------------------------------------
!> @file
!> Defines a class for handling parallelization via MPI for distributed memory architectures.
!> @author
!> Peter Bosler, Sandia National Laboratories, Albuquerque, NM
!
!> @defgroup MPISetup MPISetup
!> @brief A class for handling parallelization via MPI for distributed memory architectures.
!> 
!> The current implementation uses the 'replicated data algorithm' and simply divides the number of items by the number of processors and distributes indices accordingly.
!> Every process knows the indices assigned to every other process, and the message size associated with every other process.
!> 
!>  Note: In accordance with OpenMPI convention, these indices are assigned using a base-0 system, rather than the base-1 system typical of Fortran.
!>
!> @{
!
!
!------------------------------------------------------------------------------
use NumberKindsModule

implicit none

private
public MPISetup, New, Delete, Copy
public LoadBalance

!
!----------------
! Types and module constants
!----------------
!
type MPISetup
	integer(kint), allocatable :: indexStart(:)
	integer(kint), allocatable :: indexEnd(:)
	integer(kint), allocatable :: messageLength(:)
	integer(kint) :: n = 0
	contains
		final :: deletePrivate
end type

!
!----------------
! Interfaces
!----------------
!
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

!
!----------------
! Public methods
!----------------
!

!> @brief Allocates memory and initializes a new MPISetup object.
!> @param self Target MPISetup object
!> @param nItems The number of items to distribute across processes
!> @param nProcs The number of processes
subroutine newPrivate( self, nItems, nProcs )
	type(MPISetup), intent(out) :: self
	integer(kint), intent(in) :: nItems
	integer(kint), intent(in) :: nProcs
	
	allocate(self%indexStart(0:nProcs-1))
	allocate(self%indexEnd(0:nProcs-1))
	allocate(self%messageLength(0:nProcs-1))
	
	call LoadBalance(self, nItems, nProcs)	
end subroutine

!> @brief Deletes and frees memory associated with an MPISetup object.
!> @param self Target MPISetup object
subroutine deletePrivate(self)
	type(MPISetup), intent(inout) :: self
	if ( allocated(self%indexStart) ) then
		deallocate(self%indexStart)
	 	deallocate(self%indexEnd)
	 	deallocate(self%messageLength)
	endif
end subroutine

!> @brief Copies an MPISetup object
!> @param self Target MPISetup object
!> @param other Source MPISetup object
subroutine copyPrivate(self, other)
	type(MPISetup), intent(inout) :: self
	type(MPISetup), intent(in) :: other
	call deletePrivate(self)
	call newPrivate(self, other%n, numProcs)
end subroutine

!> @brief Distributes nItems across nProcs.
!> @param self Target MPISetup object
!> @param nItems The number of items to distribute across processes
!> @param nProcs The number of processes
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

!
!----------------
! Private methods
!----------------
!


!> @}
end module