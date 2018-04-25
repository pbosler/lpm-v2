module PSERemeshModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PlaneGeomModule
use PlanarIncompressibleModule
use BetaPlaneMeshModule
use RefinementModule
use PSEDirectSumModule
use PlanarSWEModule

implicit none

include 'mpif.h'

private

!
!----------------
! Module types and public declarations
!----------------
!

!
!----------------
! Module interfaces
!----------------
!

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PSERemesh'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains

!
!----------------
! public methods
!----------------
!

end module