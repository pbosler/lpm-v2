module MPIReplicatedDataModule
!------------------------------------------------------------------------------
! Lagrangian Particle Method
!------------------------------------------------------------------------------
!> @file
!> Defines a class for handling parallelization via MPI for distributed memory architectures.
!> @author
!> Peter Bosler, Sandia National Laboratories, Albuquerque, NM
!
!> @defgroup MPIReplicatedData MPIReplicatedData
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
use UtilitiesModule
use LoggerModule

implicit none
private

public MPIReplicatedData

type MPIReplicatedData
    integer(kint), allocatable :: startIndex(:)
    integer(kint), allocatable :: endIndex(:)
    integer(kint), allocatable :: messageSize(:)
    integer(kint) :: nItems
    contains
        procedure :: init
        procedure :: loadBalance
        final :: deleteMPI
        procedure :: logStats
end type

contains

!> @brief Allocates memory and initializes a new MPIReplicatedData object.
!> @param self Target MPIReplicatedData object
!> @param nItems The number of items to distribute across processes
!> @param nProcs The number of processes
subroutine init(self, nItems, nProcs)
    class(MPIReplicatedData), intent(inout) :: self
    integer(kint), intent(in) :: nItems
    integer(kint), intent(in) :: nProcs
    
    allocate(self%startIndex(0:nProcs-1))
    allocate(self%endIndex(0:nProcs-1))
    allocate(self%messageSize(0:nProcs-1))
    
    call self%loadBalance(nItems, nProcs)
end subroutine

!> @brief Distributes nItems across nProcs.
!> @param self Target MPIReplicatedData object
!> @param nItems The number of items to distribute across processes
!> @param nProcs The number of processes
subroutine loadBalance(self, nItems, nProcs)
    class(MPIReplicatedData), intent(inout) :: self
    integer(kint), intent(in) :: nItems, nProcs
    !
    integer(kint) :: i, chunkSize
    
    chunkSize = nItems / nProcs
    do i=0, nProcs-1
        self%startIndex(i) = i*chunkSize + 1
        self%endIndex(i) = (i+1)*chunkSize
    enddo
    self%endIndex(nProcs-1) = nItems
    self%nItems = nItems
    self%messageSize = self%endIndex - self%startIndex + 1
end subroutine

!> @brief Deletes and frees memory associated with an MPIReplicatedData object.
!> @param self Target MPIReplicatedData object
subroutine deleteMPI(self)
    type(MPIReplicatedData), intent(inout) :: self
    
    if (allocated(self%startIndex)) deallocate(self%startIndex)
    if (allocated(self%endIndex)) deallocate(self%endIndex)
    if (allocated(self%messageSize)) deallocate(self%messageSize)
end subroutine

subroutine logStats(self, aLog)
    class(MPIReplicatedData), intent(in) :: self
    type(Logger), intent(inout) :: aLog
    !
    integer(kint) :: i
    character(len=MAX_STRING_LENGTH) :: logString
    
    call StartSection(aLog, "MPIReplicatedData stats:")
    do i=1, numProcs
        write(logstring,'(A,I6, 3(A,I12))') "process rank ", i, " has indices ", self%startIndex(i), &
            " through ", self%endIndex(i), ', and messageSize = ', self%messageSize(i)
        call LogMessage(alog, TRACE_LOGGING_LEVEL, "", logstring)
    enddo
    call EndSection(alOg)
end subroutine

!>@}
end module