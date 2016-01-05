module OutputWriterModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup OutputWriter Console output writer
!> Provides an object-oriented interface for writing output to the console
!
!
! DESCRIPTION:
!> @file
!> Provides an object-oriented interface for writing output to the console
!
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
interface New
	module procedure NewPrivate
end interface

interface Delete
	module procedure DeletePrivate
end interface

interface StartSection
	module procedure StartSectionWriter
	module procedure StartBlankSectionWriter
end interface

interface EndSection
	module procedure EndSectionWriter
end interface

interface Write
	module procedure WriteString
	module procedure WriteInteger
	module procedure WriteReal
	module procedure WriteIntVector
end interface

interface WriteToMatlab
	module procedure WriteArrayToMatlab
	module procedure WriteMatrixToMatlab
	module procedure WriteIntegerMatrixToMatlab
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

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


subroutine WriteInteger(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	integer(kint), intent(in) :: val
	character(len=32) :: form
	form = FormatWithIndent(self,'(A,2X,I8)')
	write(self%fileUnit,form), trim(key),val
end subroutine

subroutine WriteReal(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: val
	character(len=32) :: form
	if ( abs(val) >= 1.0d7 ) then
		form = FormatWithIndent(self,'(A,2X,E15.8)')
	else
		form = FormatWithIndent(self,'(A,2X,F15.8)')
	endif
	write(self%fileUnit,form) trim(key), val
end subroutine

subroutine WriteIntVector(self, key, val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	type(STDIntVector), intent(in) :: val
	character(len=32) :: form
	
	write(form,'(A,I4,A)') "(A,2X,", val%N, "(I8,2X) )"
	form = FormatWithIndent(self, form)
	write(self%fileUnit, form) trim(key), val%integers(1:val%N)
end subroutine

subroutine StartBlankSectionWriter(self)
	type(OutputWriter), intent(inout) :: self
	self%indentLevel = self%indentLevel + 1
end subroutine

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


subroutine EndSectionWriter(self)
	type(OutputWriter), intent(inout) :: self
	if ( self%indentLevel == 0 ) then
		print *, "EndSection WARNING : Indentation level already 0."
	else
		self%indentLevel = self%indentLevel - 1
	endif
end subroutine

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

subroutine WriteVTKPointDataSectionHeader(fileunit, nPoints)
	integer(kint), intent(in) :: fileunit
	integer(kint), intent(in) :: nPoints
	write(fileunit,'(A,I8)') "POINT_DATA ", nPoints
end subroutine



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


!
!----------------
! Module methods : type-specific functions
!----------------
!

function FormatWithIndent(self,formatString)
	! Returns a string suitable for formatting write statements, with an indent level appropriate to
	! current state of OutputWriter::self
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: formatString
	character(len=len(formatString)+10) :: FormatWithIndent
	if ( self%indentLevel>0) then
		write(FormatWithIndent,'(A,I2,2A)') '(',TAB_SPACE*self%indentLevel,'X,',formatString(2:)
	else
		FormatWithIndent = formatString
	endif
end function


end module
