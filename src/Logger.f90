module LoggerModule
!------------------------------------------------------------------------------
! Lagrangian Particle Method
!------------------------------------------------------------------------------
!> @file
!> A Logger object for writing output to console or to files
!> @author
!> Peter Bosler, Sandia National Laboratories, Albuquerque, NM
!
!> @defgroup Logger Logger module
!> A Logger object for writing output to console or to files. 
!> 
!> Handles display, formatting, and organization of console messages.
!> All messages are presumed to have a key - value pair.  
!> For example, the key may be the origin of the message in the code, and the value contains the message content.
!> @todo Move initLogger routines from each module to this module.
!>
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

integer(kint), parameter :: DEBUG_LOGGING_LEVEL = 1 !< Priority level for debugging-related messages
integer(kint), parameter :: TRACE_LOGGING_LEVEL = 2 !< Priority level for trace-related messages and general remarks
integer(kint), parameter :: WARNING_LOGGING_LEVEL = 3 !< Priority level for warning messages
integer(kint), parameter :: ERROR_LOGGING_LEVEL = 4 !< Priority level for error messages

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
	module procedure LogIntArray
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
!> \todo In later versions, filename and fileunit will allow logging to a log file, rather than stdout
!>
!> @param self Object to be initialized
!> @param level Base priority level of messages to be logged.  Messages with lower priority than this level will be ignored.
!> @param fileunit placeholder
!> @param filename placeholder
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


!> @brief Deletes and frees memory associated with a Logger object
!> @param self Target to be deleted
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
!> @param self Target Logger object
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
!> @param self Target Logger object
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

!> @brief Logs a string/STDIntVector object
!> @param self Target Logger object
!> @param msgLevel Priority level of this message
!> @param key identification key for this message
!> @param val integer vector, message content
subroutine LogIntVector(self, msgLevel, key, val )
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	type(STDIntVector), intent(in) :: val
	
	if ( msgLevel >= self%level ) then
		call Write(self%writer, key, val)
	endif
end subroutine

subroutine LogIntArray(self, msgLevel, key, val )
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: Key
	integer(kint), dimension(:), intent(in) :: val
	
	if ( msgLevel >= self%level ) then
		call Write(self%writer, key, val )
	endif
end subroutine

!> @brief Logs a string/real key-value pair.
!> @param self Target Logger object
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
!> @param self Target Logger object
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
!> @param self Target Logger object
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
