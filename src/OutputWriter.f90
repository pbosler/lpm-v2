module OutputWriterModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
! DESCRIPTION:
!> @file
!> Provides an object-oriented interface for writing output to the console
!
!> @defgroup OutputWriter OutputWriter
!> Provides an interface for writing output to the console
!> @{
!------------------------------------------------------------------------------
use NumberKindsModule
use STDIntVectorModule

implicit none
private

public OutputWriter
public New, Delete
public StartSection, EndSection
public Write
public WriteToMatlab
public WriteVTKFileHeader, WriteVTKPointDataSectionHeader

!
!----------------
! Types and module constants
!----------------
!

type OutputWriter
	private
	integer(kint) :: fileUnit		! defines unit where output will be written
	integer(kint) :: indentLevel	! defines number of indentations to apply before writing first character
	character(len=MAX_STRING_LENGTH) :: fileName
end type

integer(KINT), parameter :: TAB_SPACE = 4

!
!----------------
! Interfaces
!----------------
!

!> @brief Initializes and OutputWriter object
interface New
	module procedure NewPrivate
end interface

!> @brief Deletes and frees memory associated with an OutputWriter object
interface Delete
	module procedure DeletePrivate
end interface

!> @brief Starts an indented section
interface StartSection
	module procedure StartSectionWriter
	module procedure StartBlankSectionWriter
end interface

!> @brief Ends an indented section
interface EndSection
	module procedure EndSectionWriter
end interface

!> @brief Writes various output to console
interface Write
	module procedure WriteString
	module procedure WriteInteger
	module procedure WriteReal
	module procedure WriteReal3
	module procedure WriteIntVector
	module procedure WriteIntegerArray
end interface

!> @brief Writes various output to a .m file for later loading in Matlab
interface WriteToMatlab
	module procedure WriteArrayToMatlab
	module procedure WriteMatrixToMatlab
	module procedure WriteIntegerMatrixToMatlab
	module procedure WriteRealScalarToMatlab
	module procedure WriteIntegerToMatlab
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

!> @brief Initializes a new OutputWriter object.
!> @todo In later versions, fileunit and filename will allow an OutputWriter's output to be redirected to a file
!>
!> @param self Target writer to be initialized
!> @param fileunit placeholder
!> @param filename placeholder
subroutine NewPrivate(self,fileUnit, fileName)
	type(OutputWriter), intent(out) :: self
	integer(kint), intent(in), optional :: fileUnit
	character(len=*), intent(in), optional :: fileName
	integer(kint) :: openstat

	self%indentLevel = 0_KINT

	if (present(fileUnit)) then
		if ( fileUnit /= STD_OUT) then
			if ( .NOT. present(fileName)) then
				self%fileUnit = STD_OUT
			else
				self%fileName = fileName
				self%fileUnit = fileUnit
			endif
		endif
	else
		self%fileUnit = STD_OUT
	endif

	if ( self%fileUnit /= STD_OUT .AND. self%fileUnit /= STD_ERR ) then
		open(unit=self%fileUnit,file=self%fileName,status='REPLACE',action='WRITE',iostat=openStat)
			if (openstat /= 0 ) then
				print *, "OutputWriter ERROR opening output file."
				print *, "Redirecting output to STD_OUT."
				self%fileUnit = STD_OUT
			endif
	endif
end subroutine

!> @brief Deletes and frees memory associated with an OutputWriter object
!> @param self Target OutputWriter object
subroutine DeletePrivate(self)
	type(OutputWriter), intent(inout) :: self
	if ( self%fileUnit /= STD_OUT .AND. self%fileUnit /= STD_ERR ) then
		close(self%fileUnit)
	endif
end subroutine

!
!----------------
! Public functions
!----------------
!

!> @brief Writes a key/value pair with appropriate indentation
!> @param self Target OutputWriter object
!> @param key identification key
!> @param str value associated with key
subroutine WriteString(self,key,str)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	character(len=*), intent(in), optional :: str
	character(len=32) :: form
	if (present(str)) then
		form = FormatWithIndent(self,'(A,2X,A)')
		write(self%fileUnit,form) trim(key), trim(str)
	else
		form = FormatWithIndent(self,'(A)')
		write(self%fileUnit,form) trim(key)
	endif
end subroutine

!> @brief Writes a key/value pair with appropriate indentation
!> @param self Target OutputWriter object
!> @param key identification key
!> @param val value associated with key
subroutine WriteInteger(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	integer(kint), intent(in) :: val
	character(len=32) :: form
	form = FormatWithIndent(self,'(A,2X,I8)')
	write(self%fileUnit,form) trim(key),val
end subroutine

!> @brief Writes a key/value pair with appropriate indentation
!> @param self Target OutputWriter object
!> @param key identification key
!> @param val value associated with key
subroutine WriteReal(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: val
	character(len=32) :: form
	if ( abs(val) >= 1.0d7 ) then
		form = FormatWithIndent(self,'(A,2X,G15.8)')
	else
		form = FormatWithIndent(self,'(A,2X,G15.8)')
	endif
	write(self%fileUnit,form) trim(key), val
end subroutine

subroutine WriteReal3(self, key, val)
    type(OutputWriter), intent(in) :: self
    character(len=*), intent(in) :: key
    real(kreal), dimension(3), intent(in) :: val
    character(len=56) :: form

    form = FormatWithIndent(self, '(A, 2X, "(", 3(G15.8,X), ")")')

    write(self%fileunit, form) trim(key), val
end subroutine

!> @brief Writes a key/value pair with appropriate indentation
!> @param self Target OutputWriter object
!> @param key identification key
!> @param val value associated with key
subroutine WriteIntVector(self, key, val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	type(STDIntVector), intent(in) :: val
	character(len=32) :: form

	write(form,'(A,I4,A)') "(A,2X,", val%N, "(I8,2X) )"
	form = FormatWithIndent(self, form)
	write(self%fileUnit, form) trim(key), val%integers(1:val%N)
end subroutine

subroutine WriteIntegerArray( self, key, val )
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	integer(kint), dimension(:), intent(in) :: val
	character(len=32) :: form
	write(form,'(A,I4,A)') "(A,2X,", size(val), "(I8,2X) )"
	form = FormatWithIndent(self,form)
	write(self%fileUnit, form) trim(key), val
end subroutine

!> @brief Starts an indented section with no title
!> @param self Target OutputWriter object
subroutine StartBlankSectionWriter(self)
	type(OutputWriter), intent(inout) :: self
	self%indentLevel = self%indentLevel + 1
end subroutine

!> @brief Starts an indented section with a title and optional description
!> @param self Target OutputWriter object
!> @param sectionName Title or subtitle of indented section
!> @param description Description of indented section
subroutine StartSectionWriter(self,sectionName,description)
	type(OutputWriter), intent(inout) :: self
	character(len=*), intent(in) :: sectionName
	character(len=*), intent(in), optional :: description
	character(len=32) :: form
	form = FormatWithIndent(self,'(A)')
	if ( procRank == 0) write(self%fileUnit,form) sectionName
	self%indentLevel = self%indentLevel + 1
	if ( present(description) ) then
		form = FormatWithIndent(self,'(A)')
		if ( procRank == 0 ) write(self%fileUnit,form) description
	endif
end subroutine

!> @brief Ends an indented sectionName
!> @param self Target OutputWriter object
subroutine EndSectionWriter(self)
	type(OutputWriter), intent(inout) :: self
	if ( self%indentLevel == 0 ) then
		print *, "EndSection WARNING : Indentation level already 0."
	else
		self%indentLevel = self%indentLevel - 1
	endif
end subroutine

!> @brief Writes the header of a Legacy formatted .vtk ASCII file
!> @param fileunit Integer unit of output .vtk file
!> @param title for .vtk file
subroutine WriteVTKFileHeader( fileunit, title )
	integer(kint), intent(in) :: fileUnit
	character(len=*), intent(in), optional :: title
	write(fileunit,'(A)') "# vtk DataFile Version 2.0"
	if ( present(title) ) then
		write(fileunit,'(A)') trim(title)
	else
		write(fileunit,'(A)') " "
	endif
	write(fileunit,'(A)') "ASCII"
	write(fileunit,'(A)') "DATASET POLYDATA"
end subroutine

!> @brief Writes the section header for VTK Point Data in a legacy formatted .vtk file
!> @param[in] fileunit Integer unit of output .vtk file
!> @param[in] nPoints Integer number of points
subroutine WriteVTKPointDataSectionHeader(fileunit, nPoints)
	integer(kint), intent(in) :: fileunit
	integer(kint), intent(in) :: nPoints
	write(fileunit,'(A,I8)') "POINT_DATA ", nPoints
end subroutine

!> @brief Writes a vector (1-d array) of real numbers to a .m file for later reading by Matlab
!> @param array Array to be written to a file
!> @param fileunit Integer unit of output .m file
!> @param name Name of array
subroutine WriteArrayToMatlab( array, fileunit, name)
	real(kreal), intent(in) :: array(:)
	integer(kint), intent(in) :: fileunit
	character(len=*), intent(in) :: name
	!
	integer(kint) :: i, n

	n = size(array)
	write(fileunit,'(A,A)',advance='NO') trim(name), " = [ "
	do i = 1, n-1
		write(fileunit,'(F24.12,A)', advance='NO') array(i), ", "
	enddo
	write(fileunit,'(F24.12,A)') array(n), "];"
end subroutine

!> @brief Writes a matrix (2-d array) of real numbers to a .m file for later reading by Matlab
!> @param matrix Matrix to be written to a file
!> @param fileunit Integer unit of output .m file
!> @param name Name of array
subroutine WriteMatrixToMatlab( matrix, fileunit, name )
	real(kreal), intent(in) :: matrix(:,:)
	integer(kint), intent(in) :: fileunit
	character(len=*), intent(in) :: name
	!
	integer(kint) :: i, j, m, n

	m = size(matrix,1)
	n = size(matrix,2)
	write(fileunit,'(A,A)',advance='NO') trim(name), " = [ "
	do i = 1, m - 1
		do j = 1, n -1
			write(fileunit, '(F24.12,A)', advance='NO') matrix(i,j), ", "
		enddo
		write(fileunit, '(F24.12,A)') matrix(i,n), "; ... "
	enddo
	do j = 1, n -1
		write(fileunit, '(F24.12,A)', advance='NO') matrix(m,j), ", "
	enddo
	write(fileunit, '(F24.12,A)') matrix(m,n), "]; "
end subroutine

!> @brief Writes a matrix (2-d array) of integers to a .m file for later reading by Matlab
!> @param matrix Matrix to be written to a file
!> @param fileunit Integer unit of output .m file
!> @param name Name of array
subroutine WriteIntegerMatrixToMatlab( matrix, fileunit, name )
	integer(kint), intent(in), dimension(:,:) :: matrix
	integer(kint), intent(in) :: fileUnit
	character(len=*), intent(in) :: name
	!
	integer(kint) :: i, j, m, n

	m = size(matrix,1)
	n = size(matrix,2)
	write(fileunit,'(A,A)', advance='NO') trim(name), " = ["
	do i = 1, m - 1
		do j = 1, n - 1
			write(fileunit,'(I8,A)', advance='no') matrix(i,j), ", "
		enddo
		write(fileunit, '(I8,A)' ) matrix(i,n), "; ..."
	enddo
	do j = 1, n - 1
		write(fileunit,'(I8,A)', advance='no') matrix(m,j), ", "
	enddo
	write(fileunit,'(I8,A)') matrix(m,n), "];"
end subroutine

subroutine WriteRealScalarToMatlab( scalar, fileUnit, name )
	real(kreal), intent(in) :: scalar
	integer(kint), intent(in) :: fileUnit
	character(len=*), intent(in) :: name

	if ( abs(scalar) >= 1.0d7 ) then
		write(fileUnit,'(A,A,E30.16,A)') trim(name), " = ", scalar, ";"
	else
		write(fileUnit,'(A,A,F30.16,A)') trim(name), " = ", scalar, ";"
	endif
end subroutine

subroutine WriteIntegerToMatlab( int, fileUnit, name )
  integer(kint), intent(in) :: int
  integer(kint), intent(in) :: fileUnit
  character(len=*), intent(in) :: name
  write(fileUnit, '(A,A,I8,A)') trim(name), " = ", int, ";"
end subroutine

!
!----------------
! Module methods : type-specific functions
!----------------
!

!> @brief This function returns a string with white space prepended to the intput formatString to provide appropriate indentation.
!> @param self Target OutputWriter object
!> @param formatString Format intsructions for console output
function FormatWithIndent(self,formatString)

	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: formatString
	character(len=len(formatString)+10) :: FormatWithIndent
	if ( self%indentLevel>0) then
		write(FormatWithIndent,'(A,I2,2A)') '(',TAB_SPACE*self%indentLevel,'X,',formatString(2:)
	else
		FormatWithIndent = formatString
	endif
end function

!> @}
end module
