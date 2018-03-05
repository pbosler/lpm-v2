module FieldOOModule

use NumberKindsModule
use LoggerModule
use ParticlesOOModule

implicit none
private
public Field

type Field
    real(kreal), allocatable :: comp1(:)
    real(kreal), allocatable :: comp2(:)
    real(kreal), allocatable :: comp3(:)
    character(MAX_STRING_LENGTH) :: name = "null"
    character(MAX_STRING_LENGTH) :: units = "null"
    integer(kint) :: N = 0
    integer(kint) :: N_Max = 0
    integer(kint) :: nDim = 0

    contains
        procedure :: init
        final :: deleteField
        procedure :: scalarMultiply
        procedure :: setToConst
        procedure :: linearComb
        procedure :: maxMagnitude
        procedure :: minMagnitude
        procedure :: maxValue
        procedure :: minValue
        procedure :: avgValue
        procedure :: logStats
        procedure :: writeMatlab
end type

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'FieldLog'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

subroutine init(self, nDim, nMax, name, units)
    class(Field), intent(inout) :: self
    integer(kint), intent(in) :: nDim, nMax
    character(len=*), intent(in), optional :: name, units

    if (.not. logInit) call InitLogger(log, procRank)

    self%N_Max = nMax
    self%N = 0
    self%nDim = nDim
    select case (nDim)
        case (1)
            allocate(self%comp1(nMax))
            self%comp1 = dzero
        case (2)
            allocate(self%comp1(nMax))
            allocate(self%comp2(nMax))
            self%comp1 = dzero
            self%comp2 = dzero
        case (3)
            allocate(self%comp1(nMax))
            allocate(self%comp2(nMax))
            allocate(self%comp3(nMax))
            self%comp1 = dzero
            self%comp2 = dzero
            self%comp3 = dzero
        case default
            call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" field init error: ", "invalid nDim.")
            return
    end select
end subroutine

subroutine deleteField(self)
    type(Field), intent(inout) :: self
    if (allocated(self%comp1)) deallocate(self%comp1)
    if (allocated(self%comp2)) deallocate(self%comp2)
    if (allocated(self%comp3)) deallocate(self%comp3)
end subroutine

pure subroutine scalarMultiply(self, a)
    class(Field), intent(inout) :: self
    real(kreal), intent(in) :: a

    select case(self%nDIm)
        case (1)
            self%comp1(1:self%N) = a*self%comp1(1:self%N)
        case (2)
            self%comp1(1:self%N) = a*self%comp1(1:self%N)
            self%comp2(1:self%N) = a*self%comp2(1:self%N)
        case (3)
            self%comp1(1:self%N) = a*self%comp1(1:self%N)
            self%comp2(1:self%N) = a*self%comp2(1:self%N)
            self%comp3(1:self%N) = a*self%comp3(1:self%N)
    end select
end subroutine

pure subroutine setToConst(self, a)
    class(Field), intent(inout) :: self
    real(kreal), intent(in) :: a
    select case (self%nDim)
        case (1)
            self%comp1(1:self%N) = a
        case (2)
            self%comp1(1:self%N) = a
            self%comp2(1:self%N) = a
        case (3)
            self%comp1(1:self%N) = a
            self%comp2(1:self%N) = a
            self%comp3(1:self%N) = a
    end select
end subroutine

subroutine linearComb(outField, a, inField1, b, infield2)
    class(Field), intent(inout) :: outField
    real(kreal), intent(in) :: a
    class(Field), intent(in) :: inField1
    real(kreal), intent(in), optional :: b
    class(Field), intent(in), optional :: inField2

    if (present(b) .and. present(infield2)) then
        select case(outField%nDim)
            case (1)
                outField%comp1(1:outField%N) = a*inField1%comp1(1:outField%N) + b*infield2%comp1(1:outfield%N)
            case (2)
                outField%comp1(1:outField%N) = a*inField1%comp1(1:outField%N) + b*infield2%comp1(1:outfield%N)
                outField%comp1(1:outField%N) = a*inField1%comp2(1:outField%N) + b*infield2%comp2(1:outfield%N)
            case (3)
                outField%comp1(1:outField%N) = a*inField1%comp1(1:outField%N) + b*infield2%comp1(1:outfield%N)
                outField%comp2(1:outField%N) = a*inField1%comp2(1:outField%N) + b*infield2%comp2(1:outfield%N)
                outField%comp3(1:outField%N) = a*inField1%comp3(1:outField%N) + b*infield2%comp3(1:outfield%N)
        end select
    else
        select case(outField%nDim)
            case (1)
                outField%comp1(1:outField%N) = a*inField1%comp1(1:outField%N)
            case (2)
                outField%comp1(1:outField%N) = a*inField1%comp1(1:outField%N)
                outField%comp1(1:outField%N) = a*inField1%comp2(1:outField%N)
            case (3)
                outField%comp1(1:outField%N) = a*inField1%comp1(1:outField%N)
                outField%comp2(1:outField%N) = a*inField1%comp2(1:outField%N)
                outField%comp3(1:outField%N) = a*inField1%comp3(1:outField%N)
        end select
    endif
end subroutine

subroutine logStats(self, aLog)
    class(Field), intent(in) :: self
    type(Logger), intent(inout) :: aLog

    call LogMessage(alog, TRACE_LOGGING_LEVEL, trim(logkey), " Field stats:")
    call StartSection(alog)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "N = ", self%N)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "N_Max = ", self%N_Max)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "minMagnitude = ", self%minMagnitude())
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "maxMagnitude = ", self%maxMagnitude())
    call EndSection(alog)
end subroutine

pure function minMagnitude(self)
    real(kreal) :: minMagnitude
    class(Field), intent(in) :: self
    !
    integer(kint) :: i
    real(kreal) :: testMag

    minMagnitude = 999.0d20
    select case (self%nDim)
        case (1)
            minMagnitude = minval(abs(self%comp1(1:self%N)))
        case (2)
            do i=1, self%N
                testMag = sqrt(self%comp1(i)**2 +self%comp2(i)**2)
                if (testMag < minMagnitude) minMagnitude = testMag
            enddo
        case (3)
            do i=1, self%N
                testMag = sqrt(self%comp1(i)**2 + self%comp2(i)**2 + self%comp3(i)**2)
                if (testMag < minMagnitude) minMagnitude = testMag
            enddo
    end select
end function

pure function maxMagnitude(self)
    real(kreal) :: maxMagnitude
    class(Field), intent(in) :: self
    integer(kint) :: i
    real(kreal) :: testMag

    maxMagnitude = dzero
    select case (self%nDim)
        case (1)
            maxMagnitude = maxval(abs(self%comp1(1:self%N)))
        case (2)
            do i=1, self%N
                testMag = sqrt(self%comp1(i)**2 +self%comp2(i)**2)
                if (testMag > maxMagnitude) maxMagnitude = testMag
            enddo
        case (3)
            do i=1, self%N
                testMag = sqrt(self%comp1(i)**2 + self%comp2(i)**2 + self%comp3(i)**2)
                if (testMag > maxMagnitude) maxMagnitude = testMag
            enddo
    end select
end function

subroutine writeMatlab(self, fileunit)
    class(Field), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i, j
    character(len=MAX_STRING_LENGTH) :: ns
    character(len=MAX_STRING_LENGTH) :: us

    if (trim(self%name) /= "null") then
        ns = trim(self%name)
    else
        ns = "field"
    endif
    write(WRITE_UNIT_1, '(A)', advance='no') trim(ns)//" = ["
    select case (self%nDim)
        case (1)
            do i=1,self%N-1
                write(WRITE_UNIT_1, '(F24.12,A)') self%comp1(i), ", "
            enddo
            write(WRITE_UNIT_1,'(F24.12,A)') self%comp1(self%N), "];"
        case (2)
            do i=1,self%N-1
                write(WRITE_UNIT_1, '(2(F24.12,A))') self%comp1(i), ", ", self%comp2(i), ";"
            enddo
            write(WRITE_UNIT_1,'(2(F24.12,A))') self%comp1(self%N), ", ", self%comp2(self%N), "];"
        case (3)
            do i=1,self%N-1
                write(WRITE_UNIT_1, '(3(F24.12,A))') self%comp1(i), ", ", self%comp2(i), ", ", self%comp3(i), ";"
            enddo
            write(WRITE_UNIT_1,'(3(F24.12,A))') self%comp1(self%N), ", ", self%comp2(self%N),", ", self%comp3(self%N), "];"
    end select
end subroutine

pure function minValue(self, component)
    real(kreal) :: minValue
    class(Field), intent(in) :: self
    integer(kint), intent(in) :: component
    minValue = 999.0d20
    select case (component)
        case (1)
            minValue = minval(self%comp1(1:self%N))
        case (2)
            minValue = minval(self%comp2(1:self%N))
        case (3)
            minValue = minval(self%comp3(1:self%N))
    end select
end function

pure function maxValue(self, component)
    real(kreal) :: maxValue
    class(Field), intent(in) :: self
    integer(kint), intent(in) :: component
    maxValue = -999.0d20
    select case (component)
        case (1)
            maxValue = maxval(self%comp1(1:self%N))
        case (2)
            maxValue = maxval(self%comp2(1:self%N))
        case (3)
            maxValue = maxval(self%comp3(1:self%N))
    end select
end function

pure function avgValue(self, component)
    real(kreal) :: avgValue
    class(Field), intent(in) :: self
    integer(kint), intent(in) :: component
    avgValue = 999.0d20
    select case(component)
        case (1)
            avgValue = sum(self%comp1(1:self%N)) / real(self%N, kreal)
        case (2)
            avgValue = sum(self%comp2(1:self%N)) / real(self%N, kreal)
        case (3)
            avgValue = sum(self%comp3(1:self%N)) / real(self%N, kreal)
    end select
end function

subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
