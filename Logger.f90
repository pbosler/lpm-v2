module LoggerModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!> @file
!> A Logger object for writing output to console or to files
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup Logger Logger module
!> A Logger object for writing output to console or to files
!> @{
!
!
!------------------------------------------------------------------------------

use NumberKindsModule
use OutputWriterModule
use STDIntVectorModule

implicit none
private

public Logger
public New, Delete
public DEBUG_LOGGING_LEVEL, TRACE_LOGGING_LEVEL, WARNING_LOGGING_LEVEL, ERROR_LOGGING_LEVEL
public LogMessage
public StartSection, EndSection

!
!----------------
! Types and module constants
!----------------
!

integer(kint), parameter :: DEBUG_LOGGING_LEVEL = 1, &
							TRACE_LOGGING_LEVEL = 2, &
							WARNING_LOGGING_LEVEL = 3, &
							ERROR_LOGGING_LEVEL = 4

!> @class Logger
!> @brief Class handles display, formatting, and organization of console messages.
!> All messages are presumed to have a key - value pair.  
!> For example, the key may be the origin of the message in the code, and the value contains the message content.
type Logger
	integer(kint) :: level !< Messages with precedence below this level will be ignored.
	type(OutputWriter) :: writer !< formatted output writer
end type

!
!----------------
! Interfaces
!----------------
!

!> @brief Initializes a logger object
interface New
	module procedure NewPrivate
end interface

!> @brief Deletes and frees memory associated with a logger object.
interface Delete
	module procedure DeletePrivate
end interface

!> @brief Common interface for logging various data types.
interface LogMessage
	module procedure LogString
	module procedure LogInteger
	module procedure LogReal
	module procedure LogIntVector
end interface

!> @brief Starts an indented section
interface StartSection
	module procedure StartSectionLogger
end interface

!> @brief Ends an indented section
interface EndSection
	module procedure EndSectionLogger
end interface

contains

!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

!> @brief Initializes a new Logger object, defines unit for chosen output.
subroutine NewPrivate(self,level,fileUnit,fileName)
	type(Logger), intent(out) :: self
	integer(kint), intent(in) :: level
	integer(kint), intent(in), optional :: fileUnit
	character(len=*), intent(in), optional :: fileName
	self%level = level

	if ( present(fileUnit) ) then
		if ( present(fileName) ) then
			call New(self%writer,fileUnit,fileName)
		else
			print *,"Logger ERROR : must supply filename for non-STD_OUT logging."
			print *,"Redirecting output to STD_OUT"
			call New(self%writer)
		endif
	else
		call New(self%writer)
	endif
end subroutine


!> @brief Placeholder for future development.
subroutine DeletePrivate(self)
	type(Logger), intent(inout) :: self
	call Delete(self%writer)
end subroutine


!
!----------------
! Public functions
!----------------
!

!> @brief Logs a string key-value pair.
!> @param self
!> @param msgLevel Priority level of this message.  
!> If it is below the priority of the logger object self, the message will not be displayed.
!> @param key identification key for message
!> @param string message content
subroutine LogString(self,msgLevel,key,string)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	character(len=*), intent(in), optional :: string
	if ( msgLevel >= self%level ) then
		if ( present(string) ) then
			call Write(self%writer,key,string)
		else
			call Write(self%writer,key)
		endif
	endif
end subroutine

!> @brief Logs a string/integer key-value pair.
!> @param self
!> @param msgLevel Priority level of this message.  
!> If it is below the priority of the logger object self, the message will not be displayed.
!> @param key identification key for message
!> @param val message content
subroutine LogInteger(self,msgLevel,key,val)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	integer(kint), intent(in) :: val
	if ( msgLevel >= self%level) then
		call Write(self%writer,key,val)
	endif
end subroutine

subroutine LogIntVector(self, msgLevel, key, val )
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	type(STDIntVector), intent(in) :: val
	
	if ( msgLevel >= self%level ) then
		call Write(self%writer, key, val)
	endif
end subroutine

!> @brief Logs a string/real key-value pair.
!> @param self
!> @param msgLevel Priority level of this message.  
!> If it is below the priority of the logger object self, the message will not be displayed.
!> @param key identification key for message
!> @param val message content
subroutine LogReal(self,msgLevel,key,val)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: val
	if (msgLevel>=self%level) then
		call Write(self%writer,key,val)
	endif
end subroutine

!> @brief Starts an indented section for this logger object.
!> Allows user to optionally provide a section title and section description.
!> @param self
!> @param sectionName
!> @param description
subroutine StartSectionLogger(self,sectionName,description)
	type(Logger), intent(inout) :: self
	character(len=*), intent(in), optional :: sectionName
	character(len=*), intent(in), optional :: description
	
	if ( present(sectionName) ) then
		if ( present(description) ) then
			call StartSection(self%writer,sectionName,description)
		else
			call StartSection(self%writer,sectionName)
		endif
	else
		call StartSection(self%writer)
	endif
end subroutine


!> @brief Ends an indented section for this logger object.
!> @param self
subroutine EndSectionLogger(self)
	type(Logger), intent(inout) :: self
	call EndSection(self%writer)
end subroutine


!
!----------------
! Module methods : type-specific functions
!----------------
!


!> @}
end module
