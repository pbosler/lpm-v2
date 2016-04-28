module FieldModule
!------------------------------------------------------------------------------
! Lagrangian Particle Method (LPM) version 1.5
!------------------------------------------------------------------------------
!> @file
!> Provides a data type for holding variables associated with an LPM spatial discretization by instances of a Particles object.
!> @author
!> Peter Bosler, Sandia National Laboratories Center for Computing Research
!
!> @defgroup Field Field
!> @brief Provides a vectorized data structure for scalar and vector fields defined on @ref Particles objects.
!>
!> Allows users to set a name and units, if applicable, for each field.
!> Provides a data type for holding variables associated with an LPM spatial discretization by instances of @ref Particles .
!> 
!> May be scalar or vector data.
!> Field objects have a 1-to-1 correspondence with a @ref Particles object, so that the data associated with 
!> the particle whose index is i are in the field object also at index i. 
!>
!> The field data structures is a "structure of arrays," so that all information about a scalar of vector field at
!> particle i is located at index i in the field arrays.  
!> For example, the value of a scalar carried by particle i is `aField%%scalar(i)`. 
!> Vector field components must be accessed individually using `aField%%xComp(i)`, `aField%%yComp(i)`, etc.
!>
!>
!> @{
!
!------------------------------------------------------------------------------
use NumberKindsModule
use LoggerModule
use ParticlesModule
use FacesModule

implicit none
private
public Field
public New, Delete, Copy
public InsertScalarToField, InsertVectorToField
public WriteFieldToVTKPointData, WriteFieldToVTKCellData
public WriteFieldToMatlab
public LogStats
public SetFieldToZero
public MinMagnitude, MaxMagnitude, MaxScalarVal, MinScalarVal
public MultiplyFieldByScalar
public ScalarAverage
!public SetFieldToScalarFunction, SetFieldToVectorFunction

!> @brief Vectorized class for physical data defined on a particles' spatial discretization of a domain.
!> May be scalar or vector data.
!>
!> Field objects have a 1-to-1 correspondence with a @ref Particles object, so that the data associated with 
!> the particle whose index is i are in the field object also at index i. 
type Field
	real(kreal), allocatable :: scalar(:) !< Array to hold scalar field data
	real(kreal), allocatable :: xComp(:)  !< Array to hold the x-component of vector field data
	real(kreal), allocatable :: yComp(:)  !< Array to hold the y-component of vector field data
	real(kreal), allocatable :: zComp(:)  !< Array to hold the z-component of vector field data
	character(MAX_STRING_LENGTH) :: name = "null" !< name of field variable
	character(MAX_STRING_LENGTH) :: units = "null" !< physical units of field variable
	integer(kint) :: N = 0 !< Number of field values stored (should be equal to the number of particles)
	integer(kint) :: N_Max = 0 !< maximum number of field values allowed in memory
	integer(kint) :: nDim = 0 !< dimension of field, 1 = scalar, >= 2 for vector fields
	
	contains	
		final :: deletePrivate
end type	

!> @brief Allocates memory and initializes to null/zero a Field object.
interface New
	module procedure newPrivate
end interface

!> @brief Deletes a Field object and frees its memory.
interface Delete
	module procedure deletePrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

!> @brief Outputs statistics about a Field object to the console via a @ref Logger object.
interface LogStats
	module procedure logStatsPrivate
end interface
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'FieldLog'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

!> @brief Allocates memory and initializes to zero a Field object.
!> @param[out] self Target field object
!> @param[in] nDim number of components in this Field object (not the spatial domain of its accompanying particle set)
!> @param[in] nMax max number of particles that can be associated with this field
!> @param[in] name e.g., vorticity or potential
!> @param[in] units physical units, if applicable
subroutine newPrivate(self, nDim, nMax, name, units )
	type(Field), intent(out) :: self
	integer(kint), intent(in) :: nDim, nMax
	character(len=*), intent(in), optional :: name, units
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	if ( nDim <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " invalid nMax.")
		return
	endif
	
	self%N_Max = nMax
	self%nDim = nDim
	self%N = 0
	select case (nDim)
		case (1)
			allocate(self%scalar(nMax))
			self%scalar = 0.0
		case (2)
			allocate(self%xComp(nMax))
			allocate(self%yComp(nMax))
			self%xComp = 0.0
			self%yComp = 0.0
		case (3)
			allocate(self%xComp(nMax))
			allocate(self%yComp(nMax))
			allocate(self%zComp(nMax))
			self%xComp = 0.0
			self%yComp = 0.0
			self%zComp = 0.0
		case default
	end select 
	
	if ( present(name) ) then
		self%name = trim(name)
	endif
	
	if ( present(units) ) then
		self%units = trim(units)
	endif
end subroutine

!> @brief Deallocates memory assigned by newprivate.
!> @param self Target field object.
subroutine deletePrivate(self)
	type(Field), intent(inout) :: self
	if ( allocated(self%scalar) ) deallocate(self%scalar)
	if ( allocated(self%xComp) ) deallocate(self%xComp)
	if ( allocated(self%yComp) ) deallocate(self%yComp)
	if ( allocated(self%zComp) ) deallocate(self%zComp)
end subroutine

!> @brief Performs a deep copy of one field object into another.
!> @param self Target field object
!> @param other Source field object
subroutine copyPrivate( self, other )
	type(Field), intent(inout) :: self
	type(Field), intent(in) :: other
	integer(kint) :: i
	if ( self%nDim /= other%nDim .OR. self%N_Max < other%N ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " Copy Field ERROR : memory not allocated.") 
		return 
	endif
	self%N = other%N
	self%name = other%name
	self%units = other%units
	if ( self%nDim == 1 ) then
		do i = 1, other%N
			self%scalar(i) = other%scalar(i)
		enddo
		self%scalar(self%N+1:self%N_Max) = 0.0_kreal
	elseif ( self%nDim == 2 ) then
		do i = 1, other%N
			self%xComp(i) = other%xComp(i)
			self%yComp(i) = other%yComp(i)
		enddo
		
	else
		do i = 1, other%N
			self%xComp(i) = other%xComp(i)
			self%yComp(i) = other%yComp(i)
			self%zComp(i) = other%zComp(i)
		enddo
	endif
end subroutine

!> @brief Inserts a scalar value to a preallocated Field object.
!> @param self Target field object
!> @param val 
subroutine InsertScalarToField(self, val)
	type(Field), intent(inout) :: self
	real(kreal), intent(in) :: val
	if ( .NOT. allocated(self%scalar) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" "//self%name," InsertScalarToField : scalar not allocated.")
		return
	endif
	if ( self%N >= self%N_Max ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" "//self%name, " InsertScalarToField : out of memory.")
		return
	endif
	self%scalar( self%N + 1 ) = val
	self%N = self%N + 1
end subroutine

!> @brief Inserts a vector value to a preallocated Field object.
!> @param self Target field object
!> @param vecval
subroutine InsertVectorToField( self, vecVal )
	type(Field), intent(inout) :: self
	real(kreal), intent(in) :: vecVal(:)
	
	if ( self%N >= self%N_Max ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertVectorToField : out of memory.")
		return
	endif
	if ( self%nDim == 2 .AND. size(vecVal) == 2 ) then
		self%xComp( self%N + 1 ) = vecVal(1)
		self%yComp( self%N + 1 ) = vecVal(2)
	elseif ( self%nDim == 3 .AND. size(vecVal) == 3 ) then
		self%xComp( self%N + 1 ) = vecVal(1)
		self%yComp( self%N + 1 ) = vecVal(2)
		self%zComp( self%N + 1 ) = vecVal(3)
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertVectorToField : size mismatch.")
		return
	endif
	self%N = self%N + 1
end subroutine


!> @brief Finds the area-average of a scalar.
!> @param self Target field object
!> @param aParticles particlesmodule::particles object associated with this field
!> @return @f$ \frac{1}{A} \int f(x)\, dA @f$
function ScalarAverage(self, aParticles )
	real(kreal) :: ScalarAverage
	type(Field), intent(in) :: self
	type(Particles), intent(in) :: aParticles
	ScalarAverage = sum( self%scalar(1:self%N) * aParticles%area(1:self%N), MASK=aParticles%isActive(1:self%N)) / &
					TotalArea(aParticles)
end function

!> @brief Zeroes all values in a field object.
!> @param self Target field object
subroutine SetFieldToZero( self )
	type(Field), intent(inout) :: self

	if (self%nDim == 1) then
		self%scalar = 0.0_kreal
	elseif ( self%nDim == 2 ) then
		self%xComp = 0.0_kreal
		self%yComp = 0.0_kreal
	elseif ( self%nDim == 3 ) then
		self%xComp = 0.0_kreal
		self%yComp = 0.0_kreal
		self%zComp = 0.0_kreal
	endif
end subroutine

!> @brief Multiplies all values in a field by a scalar
!> @param self Target field object
!> @param multiplier
subroutine MultiplyFieldByScalar(self, multiplier)
	type(Field), intent(inout) :: self
	real(kreal), intent(in) :: multiplier
	if (self%nDim == 1) then
		self%scalar = multiplier * self%scalar
	elseif ( self%nDim == 2 ) then
		self%xComp = multiplier * self%xComp
		self%yComp = multiplier * self%yComp
	elseif ( self%nDim == 3 ) then
		self%xComp = multiplier * self%xComp
		self%yComp = multiplier * self%yComp
		self%zComp = multiplier * self%zComp
	endif
end subroutine

!> @brief Outputs Field object data to a .vtk PolyData file for use with VTK or ParaView
!> @param self
!> @param fileunit
subroutine WriteFieldToVTKPointData( self, fileunit )
	type(Field), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	!
	integer(kint) :: j
	
	!write(fileunit,'(A,I8)') "POINT_DATA ", self%N
	write(fileunit,'(5A,I4)') "SCALARS ", trim(self%name), "_", trim(self%units) , "  double ", self%nDim
	write(fileunit,'(A)') "LOOKUP_TABLE default"
	
	select case (self%nDim)
		case (1)
			do j = 1, self%N
				if ( abs(self%scalar(j)) < ZERO_TOL ) then
					write(fileunit, *) 0.0_kreal
				else
					write(fileunit,*) self%scalar(j)
				endif
			enddo
		case (2)
			do j = 1, self%N
				write(fileunit,*) self%xComp(j), self%yComp(j)
			enddo
		case (3)
			do j = 1, self%N
				write(fileunit,*) self%xComp(j), self%yComp(j), self%zComp(j)
			enddo
	end select
end subroutine

!> @brief Outputs Field object data to a .vtk PolyData file for use with VTK or ParaView, 
!> if this field is also associated with a faces object.
!> @param self
!> @param fileunit
!> @param aFaces
subroutine WriteFieldToVTKCellData( self, fileunit, aFaces )
	type(Field), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	type(Faces), intent(in) :: aFaces
	!
	integer(kint) :: i, j, nCells, nVerts
	
	nVerts = size(aFaces%vertices, 1)
	nCells = nVerts * aFaces%N_Active
	
	!write(fileunit,'(A,I8)') "CELL_DATA ", nCells
	write(fileunit,'(5A,I4)') "SCALARS ", trim(self%name), "_", trim(self%units) , "  double ", self%nDim
	write(fileunit,'(A)') "LOOKUP_TABLE default"
	
	select case ( self%nDim )
		case (1)
			do i = 1, aFaces%N
				if ( .NOT. aFaces%hasChildren(i) ) then
					do j = 1, nVerts
						write(fileunit,*) self%scalar( aFaces%centerParticle(i) )
					enddo
				endif
			enddo 
		case (2)
			do i = 1, aFaces%N
				if ( .NOT. aFaces%hasChildren(i) ) then
					do j = 1, nVerts
						write(fileunit, *) self%xComp( aFaces%centerParticle(i)), self%yComp( aFaces%centerParticle(j))
					enddo
				endif
			enddo
		case (3)
			do i = 1, aFaces%N
				if ( .NOT. aFaces%hasChildren(i) ) then
					do j = 1, nVerts
						write(fileunit, *) self%xComp( aFaces%centerParticle(j)), self%yComp( aFaces%centerParticle(j)),&
									   self%zComp(aFaces%centerParticle(j))
					enddo
				endif
			enddo
	end select
end subroutine

!> @brief Writes Field information to console using a loggermodule::logger object for formatting. 
!> @param self
!> @param aLog
subroutine logStatsPrivate(self, aLog )
	type(Field), intent(in) :: self
	type(Logger), intent(inout) :: aLog 
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, logkey, "FieldStats:"//trim(self%name) )
	call StartSection(aLog)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "Field.nDim = ", self%nDim)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "Field.N = ", self%N)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "Field.name = ", self%name)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "Field.units = ", self%units)
	if (self%nDim == 1 ) then
		 call LogMessage(alog, TRACE_LOGGING_LEVEL,"max scalar = ", maxval(self%scalar(1:self%N)))
		 call LogMessage(alog, TRACE_LOGGING_LEVEL,"min scalar = ", minval(self%scalar(1:self%N)))
	else
		call LogMessage(alog, TRACE_LOGGING_LEVEL, "max magnitude = ", MaxMagnitude(self) ) 
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "min magnitude = ", MinMagnitude(self) )
	endif
	call EndSection(aLog)
end subroutine

!> @brief Returns the maximum magnitude of a vector field
!> @param self
pure function MaxMagnitude(self)
	real(kreal) :: MaxMagnitude 
	type(Field), intent(in) :: self
	!
	integer(kint) :: i
	real(kreal) :: mag
	MaxMagnitude = -1.0d64
	if (self%nDim == 2 ) then
		do i = 1, self%N
			mag = sqrt( self%xComp(i)*self%xComp(i) + self%yComp(i)*self%yComp(i))
			if ( mag > MaxMagnitude ) MaxMagnitude = mag
		enddo
	else if ( self%nDim == 3 ) then
		do i = 1, self%N
			mag = sqrt( self%xComp(i)*self%xComp(i) + self%yComp(i)*self%yComp(i) + self%zComp(i)*self%zComp(i))
			if ( mag > MaxMagnitude ) MaxMagnitude = mag
		enddo 	
	endif
end function 

!> @brief Returns the minimum magnitude of a vector field
!> @param self
pure function MinMagnitude(self)
	real(kreal) :: MinMagnitude 
	type(Field), intent(in) :: self
	!
	integer(kint) :: i
	real(kreal) :: mag
	MinMagnitude = 1.0d64
	if (self%nDim == 2 ) then
		do i = 1, self%N
			mag = sqrt( self%xComp(i)*self%xComp(i) + self%yComp(i)*self%yComp(i))
			if ( mag < MinMagnitude ) MinMagnitude = mag
		enddo
	else if ( self%nDim == 3 ) then
		do i = 1, self%N
			mag = sqrt( self%xComp(i)*self%xComp(i) + self%yComp(i)*self%yComp(i) + self%zComp(i)*self%zComp(i))
			if ( mag < MinMagnitude ) MinMagnitude = mag
		enddo 	
	endif
end function

pure function MaxScalarVal(self)
	real(kreal) :: MaxScalarVal
	type(Field), intent(in) :: self
	MaxScalarVal = maxval( self%scalar(1:self%N) )
end function

pure function MinScalarVal(self)
	real(kreal) :: MinScalarVal
	type(Field), intent(in) :: self
	MinScalarVal = minval( self%scalar(1:self%N) )
end function

!> @brief Writes Field data to a script .m file readable by Matlab
!> @param self
!> @param fileunit
subroutine WriteFieldToMatlab( self, fileunit )
	type(Field), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	!
	integer(kint) :: i
	
	select case ( self%nDim )
		case (1)
			write(fileunit,*) "scalarField_", trim(self%name), " = [", self%scalar(1), ", ..."
			do i = 2, self%N - 1
				write(fileunit, *) self%scalar(i), ", ..."
			enddo
			write(fileunit, *) self%scalar(self%N), "];"
		case (2)
			write(fileunit,*) "vectorField_", trim(self%name), " = [", self%xComp(1), ", ", self%yComp(1), "; ..."
			do i = 2, self%N - 1
				write(fileunit,*) self%xComp(i), ", ", self%yComp(i), "; ..."
			enddo
			write(fileunit,*) self%xComp(self%N), ", ", self%yComp(self%N), "];"
		case (3)
			write(fileunit,*) "vectorField_", trim(self%name), " = [", self%xComp(1), ", ",&
										 self%yComp(1), ", ", self%zComp(1), "; ..."
			do i = 2, self%N-1
				write(fileunit,*) self%xComp(i), ", ", self%yComp(i), ", ", self%zComp(i), "; ..."
			enddo							 
			write(fileunit,*) self%xComp(self%N), ", ", self%yComp(self%N), ", ", self%zComp(self%N), "];"
	end select	
end subroutine

!
!----------------
! Private methods
!----------------
!


!> @brief Initializes a logger for the Field module
!> 
!> Output is controlled both by message priority and by MPI Rank
!> @param aLog Target Logger object
!> @param rank Rank of this processor
subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

!> @}
end module