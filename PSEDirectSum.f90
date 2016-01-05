module PSEDirectSumModule

use NumberKindsModule
use LoggerModule
use ParticlesModule
use EdgesModule, only : MaxEdgeLength
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use SphereGeomModule, only : ChordDistance, SphereDistance, SphereProjection

implicit none

include 'mpif.h'

private
public PSE, New, Delete
public PSEPlaneInterpolateScalar
public PSEPlaneGradientAtParticles, PSEPlaneSecondPartialsAtParticles, PSEPlaneLaplacianAtParticles
public PSESphereInterpolateScalar, PSESphereGradientAtParticles, PSESphereDivergenceAtParticles, PSESphereLaplacianAtParticles
public bivariateFirstDerivativeKernel8, bivariateLaplacianKernel8
public PSEPlaneDoubleDotProductAtParticles

type PSE
	real(kreal) :: eps
	contains
		final :: deletePrivate
end type

interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PSE'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!----------------
!
! Public functions
!
!----------------
subroutine newPrivate( self, aMesh, radiusMultiplier )
	type(PSE), intent(out) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in), optional :: radiusMultiplier
	!
	real(kreal) :: multiple
	integer(kint) :: nP
	
	if ( .NOT. logInit ) call InitLogger(log, procRank) 
	
	if ( present(radiusMultiplier) ) then
		multiple = radiusMultiplier
	else
		multiple = 2.0_kreal
	endif
	
	self%eps = multiple * MaxEdgeLength( aMesh%edges, aMesh%particles)
	
end subroutine

subroutine deletePrivate(self)
	type(PSE), intent(inout) :: self
end subroutine

pure function PSEPlaneInterpolateScalar(self, aMesh, scalarField, interpLoc )
	real(kreal) :: PSEPlaneInterpolateScalar
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	real(kreal), intent(in) :: interpLoc(3)
	!
	integer(kint) :: j
	real(kreal) :: kIn, xj(3)
	
	PSEPlaneInterpolateScalar = 0.0_kreal
	do j = 1, aMesh%particles%N
		if ( aMesh%particles%isActive(j) ) then
			xj = PhysCoord(aMesh%particles,j)
			kIn = ChordDistance(interpLoc,xj)/self%eps
			PSEPlaneInterpolateScalar = PSEPlaneInterpolateScalar + scalarField%scalar(j) * aMesh%particles%area(j) * &
			 	bivariateDeltaKernel8( kIn ) / (self%eps**2)
		endif
	enddo
end function


subroutine PSEPlaneGradientAtParticles( self, aMesh, scalarField, scalarGrad, particlesMPI )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarGrad
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3), kOut
	
	call SetFieldToZero(scalarGrad)
	scalarGrad%N = aMesh%particles%N
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey), " entering PSEGradientAtParticles...")
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = ChordDistance(xi,xj)/self%eps
				scalarGrad%xComp(i) = scalarGrad%xComp(i) + ( scalarField%scalar(j) + scalarField%scalar(i)) *&
					(xi(1)-xj(1)) * bivariateFirstDerivativeKernel8(kIn) / (self%eps**3) * aMesh%particles%area(j)
				scalarGrad%yComp(i) = scalarGrad%yComp(i) + ( scalarField%scalar(j) + scalarField%scalar(i)) *&
					(xi(2)-xj(2)) * bivariateFirstDerivativeKernel8(kIn) / (self%eps**3) * aMesh%particles%area(j)					
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(scalarGrad%xComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
					   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(scalarGrad%yComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
					   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)					   
	enddo
	
	call MultiplyFieldByScalar(scalarGrad, 1.0_kreal/self%eps)
end subroutine

subroutine PSESphereGradientAtParticles( self, amesh, scalarField, scalarGrad, particlesMPI )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarGrad
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3), grad(3), proj(3,3), kOut
	
	call SetFieldToZero(scalarGrad)
	scalarGrad%N = aMesh%particles%N
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey), " entering PSESphereGradientAtParticles.")
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		proj = SphereProjection(xi)
		grad = 0.0_kreal
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord( aMesh%particles, j)
				kIn = SphereDistance(xj, xi) / self%eps
				kOut = bivariateFirstDerivativeKernel8( kIn ) / self%eps**2
				grad(1) = grad(1) + ( scalarField%scalar(j) + scalarField%scalar(i)) * ( xi(1) - xj(1) ) * &
						  kOut * aMesh%particles%area(j)
				grad(2) = grad(2) + ( scalarField%scalar(j) + scalarField%scalar(i)) * ( xi(2) - xj(2) ) * &
						  kOut * aMesh%particles%area(j)
			    grad(3) = grad(3) + ( scalarField%scalar(j) + scalarField%scalar(i)) * ( xi(3) - xj(3) ) * &
			    		  kOut * aMesh%particles%area(j)
			endif
		enddo
		grad = MATMUL( proj, grad)
		scalarGrad%xComp(i) = grad(1)
		scalarGrad%yComp(i) = grad(2)
		scalarGrad%zComp(i) = grad(3)
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( scalarGrad%xComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
			particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( scalarGrad%yComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
			particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( scalarGrad%zComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
			particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)	
	enddo
	call MultiplyFieldByScalar(scalarGrad, 1.0_kreal/(self%eps**2))
end subroutine

subroutine PSEPlaneSecondPartialsAtParticles( self, aMesh, scalarGrad, secondPartials, particlesMPI )
	type(PSE), intent(in) :: self
	class(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarGrad
	type(Field), intent(inout) :: secondPartials
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3), dxx, dxy, dyx, dyy, kOut
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey), " entering PSESecondPartialsAtParticles...")
	
	call SetFieldToZero(secondPartials)
	secondPartials%N = aMesh%particles%N
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord( aMesh%particles, i)
		dxx = 0.0_kreal
		dxy = 0.0_kreal
		dyx = 0.0_kreal
		dyy = 0.0_kreal
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = ChordDistance(xi,xj)/self%eps
				kOUt = bivariateFirstDerivativeKernel8(kIn) / self%eps**2
				dxx = dxx + (scalarGrad%xComp(j) + scalarGrad%xComp(i)) * ( xi(1) - xj(1) )/self%eps * &
					kOut * aMesh%particles%area(j)
				dxy = dxy + (scalarGrad%xComp(j) + scalarGrad%xComp(i)) * ( xi(2) - xj(2) )/self%eps * &
					kOut * aMesh%particles%area(j)
				dyx = dyx + (scalarGrad%yComp(j) + scalarGrad%yComp(i)) * ( xi(1) - xj(1)) /self%eps * &
					kOut * aMesh%particles%area(j)
				dyy = dyy + (scalarGrad%yComp(j) + scalarGrad%yComp(i)) * ( xi(2) - xj(2)) /self%eps * &
					kOut * aMesh%particles%area(j)
			endif
		enddo
		secondPartials%xComp(i) = dxx
		secondPartials%yComp(i) = 0.5_kreal * ( dxy + dyx )
		secondPartials%zComp(i) = dyy
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(secondPartials%xComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
				   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCOde)
		call MPI_BCAST(secondPartials%yComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
				   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCOde)
		call MPI_BCAST(secondPartials%zComp(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
				   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCOde)					   
	enddo
	
	call MultiplyFieldByScalar( secondPartials, 1.0_kreal / self%eps )
end subroutine

subroutine PSEPlaneDoubleDotProductAtParticles( self, aMesh, vectorField, doubleDot, particlesMPI)
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: vectorField
	type(Field), intent(inout) :: doubleDot
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3), ux, uy, vx, vy, kOut
	
	call SetFieldToZero(doubleDot)
	doubleDot%N = aMesh%particles%N
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey), " entering PSEDoubleDotProduct ...")
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		ux = 0.0_kreal
		uy = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = ChordDistance(xi,xj)/self%eps
				kOut = bivariateFirstDerivativeKernel8( kIn ) / self%eps**2
				ux = ux + (vectorField%xComp(j) + vectorField%xComp(i)) * ( xi(1) - xj(1) )/self%eps * &
					kOut * aMesh%particles%area(j)
				uy = uy + (vectorField%xComp(j) + vectorField%xComp(i)) * ( xi(2) - xj(2) )/self%eps * &
					kOut * aMesh%particles%area(j)
				vx = vx + (vectorField%yComp(j) + vectorField%yComp(i)) * ( xi(1) - xj(1) )/self%eps * &
					kOut * aMesh%particles%area(j)	
				vy = vy + (vectorField%yComp(j) + vectorField%yComp(i)) * ( xi(2) - xj(2) )/self%eps * &
					kOut * aMesh%particles%area(j)		
			endif
		enddo
		doubleDot%scalar(i) = (ux * ux + 2.0_kreal * uy * vx + vy * vy)/self%eps**2
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( doubleDot%scalar(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
				particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
	enddo	
end subroutine

subroutine PSESphereDoubleDotProductAtParticles( self, aMesh, vectorField, doubleDot, particlesMPI)
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: vectorField
	type(Field), intent(inout) :: doubleDot
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3), ux, uy, uz, vx, vy, vz, wx, wy, wz, kOut
	
	call SetFieldToZero(doubleDot)
	doubleDot%N = aMesh%particles%N
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		ux = 0.0_kreal
		uy = 0.0_kreal
		uz = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		vz = 0.0_kreal
		wx = 0.0_kreal
		wy = 0.0_kreal
		wz = 0.0_kreal
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = SphereDistance(xi,xj) / self%eps
				kOut = bivariateFirstDerivativeKernel8( kIn ) / self%eps**2
				ux = ux + (vectorField%xComp(j) + vectorField%xComp(i)) * ( xi(1) - xj(1) )/self%eps * &
					kOut * aMesh%particles%area(j)
				uy = uy + (vectorField%xComp(j) + vectorField%xComp(i)) * ( xi(2) - xj(2) )/self%eps * &
					kOut * aMesh%particles%area(j)
				uz = uz + (vectorField%xComp(j) + vectorField%xComp(i)) * ( xi(3) - xj(3) )/self%eps * &
					kOut * aMesh%particles%area(j)
				vx = vx + (vectorField%yComp(j) + vectorField%yComp(i)) * ( xi(1) - xj(1) )/self%eps * &
					kOut * aMesh%particles%area(j)
				vy = vy + (vectorField%yComp(j) + vectorField%yComp(i)) * ( xi(2) - xj(2) )/self%eps * &
					kOut * aMesh%particles%area(j)
				vz = vz + (vectorField%yComp(j) + vectorField%yComp(i)) * ( xi(3) - xj(3) )/self%eps * &
					kOut * aMesh%particles%area(j)	
				wx = wx + (vectorField%zComp(j) + vectorField%yComp(i)) * ( xi(1) - xj(1) )/self%eps * &
					kOut * aMesh%particles%area(j)
				wy = wy + (vectorField%zComp(j) + vectorField%yComp(i)) * ( xi(2) - xj(2) )/self%eps * &
					kOut * aMesh%particles%area(j)
				wz = wz + (vectorField%zComp(j) + vectorField%yComp(i)) * ( xi(3) - xj(3) )/self%eps * &
					kOut * aMesh%particles%area(j)	
			endif
		enddo
		doubleDot%scalar(i) = (ux*ux + vy*vy + wz*wz + 2.0_kreal*( uy*vx + uz*wx + vz*wy))/self%eps**2
	enddo
end subroutine

subroutine PSEDoubleDotProductArrays( self, doubleDot, xIn, yIn, uIn, vIn, areaIn, activeMask, mpiParticles )
	type(PSE), intent(in) :: self
	real(kreal), intent(out) :: doubleDot(:)
	real(kreal), intent(in) :: xIn(:)
	real(kreal), intent(in) :: yIn(:)
	real(kreal), intent(in) :: uIn(:)
	real(kreal), intent(in) :: vIn(:)
	real(kreal), intent(in) :: areaIn(:)
	logical(klog), intent(in) :: activeMask(:)
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, ux, uy, vx, vy, kernel
	
	doubleDot(mpiParticles%indexStart(procRank):mpiParticles%indexEnd(procRank)) = 0.0_kreal
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		ux = 0.0_kreal
		uy = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		if ( activeMask(j) ) then
			kIn  = sqrt( (xIn(i) - xIn(j))**2 + (yIn(i) - yIn(j))**2 ) / self%eps
			kernel = bivariateFirstDerivativeKernel8( kIn ) / self%eps**2
			ux = ux + ( uIn(j) + uIn(i) ) * ( xIn(i) - xIn(j) ) / self%eps * kernel * areaIn(j)
			uy = uy + ( uIn(j) + uIn(i) ) * ( yIn(i) - yIn(j) ) / self%eps * kernel * areaIn(j)
			vx = vx + ( vIn(j) + vIn(i) ) * ( xIn(i) - xIn(j) ) / self%eps * kernel * areaIn(j)
			vy = vy + ( vIn(j) + vIn(i) ) * ( yIn(i) - yIn(j) ) / self%eps * kernel * areaIn(j)				
		endif
		doubleDot(i) = ( ux * ux + 2.0_kreal * uy * vx + vy * vy ) / self%eps**2
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( doubleDot(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
end subroutine

subroutine PSEPlaneLaplacianAtParticles(self, aMesh, scalarField, scalarLap, particlesMPI )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarLap
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3)
	
	call SetFieldToZero(scalarLap)
	scalarLap%N = aMesh%particles%N
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey), " entering PSELaplacianAtParticles...")
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = ChordDistance(xi,xj)/self%eps
				scalarLap%scalar(i) = scalarLap%scalar(i) + (scalarField%scalar(j) - scalarField%scalar(i)) * &
					bivariateLaplacianKernel8(kIN)/(self%eps**2) * aMesh%particles%area(j)
			endif
		enddo
	enddo
	
	do i = 0, numProcs-1
		call MPI_BCAST(scalarLap%scalar(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
					   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
	
	call MultiplyFieldByScalar(scalarLap, 1.0_kreal/(self%eps**2))
end subroutine

subroutine PSESphereLaplacianAtParticles( self, aMesh, scalarField, scalarLap, particlesMPI )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarLap
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3)
	
	call SetFieldToZero(scalarLap)
	scalarLap%N = aMesh%particles%N
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey), " entering PSELaplacianAtParticles...")
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = SphereDistance(xi,xj)/self%eps
				scalarLap%scalar(i) = scalarLap%scalar(i) + (scalarField%scalar(j) - scalarField%scalar(i)) * &
					bivariateLaplacianKernel8(kIN)/(self%eps**2) * aMesh%particles%area(j)
			endif
		enddo
	enddo
	
	do i = 0, numProcs-1
		call MPI_BCAST(scalarLap%scalar(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
					   particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
	
	call MultiplyFieldByScalar(scalarLap, 1.0_kreal/(self%eps**2))
end subroutine

pure function PSESphereInterpolateScalar( self, aMesh, scalarField, interpLoc)
	real(kreal) :: PSESphereInterpolateScalar
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	real(kreal), intent(in) :: interpLoc(3)
	!
	integer(kint) :: j
	real(kreal) :: kIn, xj(3)
	
	PSESphereInterpolateScalar = 0.0_kreal
	do j = 1, aMesh%particles%N
		if ( aMesh%particles%isActive(j) ) then
			xj = PhysCoord(aMesh%particles, j)
			kIn = SphereDistance(xj, interpLoc)/self%eps
			PSESphereInterpolateScalar = PSESphereInterpolateScalar + scalarField%scalar(j) * aMesh%particles%area(j) * &
				bivariateDeltaKernel8( kIn ) / (self%eps**2)
		endif
	enddo
end function

subroutine PSESphereDivergenceAtParticles( self, aMesh, vectorField, divergence, particlesMPI )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: vectorField
	type(Field), intent(inout) :: divergence
	type(MPISetup), intent(in) :: particlesMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, xi(3), xj(3), grad(3), proj(3,3), denom
	
	call SetFieldToZero(divergence)
	divergence%N = aMesh%particles%N
	
	denom = self%eps**3
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey), " entering PSESphereDivergenceAtParticles.")
	
	do i = particlesMPI%indexStart(procRank), particlesMPI%indexEnd(procRank)
		xi = PhysCoord(aMesh%particles, i)
		proj = SphereProjection(xi)
		grad = 0.0_kreal
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				kIn = SphereDistance(xj, xi)/self%eps
				grad(1) = (xi(1) - xj(1)) * bivariateFirstDerivativeKernel8( kIn ) / denom
				grad(2) = (xi(2) - xj(2)) * bivariateFirstDerivativeKernel8( kIn ) / denom
				grad(3) = (xi(3) - xj(3)) * bivariateFirstDerivativeKernel8( kIn ) / denom
				grad = MATMUL( proj, grad )
				divergence%scalar(i) = divergence%scalar(i) + ( grad(1) * (vectorField%xComp(j) + vectorField%xComp(i)) + &
						grad(2) * (vectorField%yComp(j) + vectorField%yComp(i) ) +&
						grad(3) * (vectorField%zComp(j) + vectorField%zComp(i) )) * aMesh%particles%area(j)
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( divergence%scalar(particlesMPI%indexStart(i):particlesMPI%indexEnd(i)), &
			particlesMPI%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
	enddo
	
	call MultiplyFieldByScalar(divergence, 1.0_kreal/self%eps)
end subroutine

!----------------
!
! Private functions
!
!----------------

pure function integralKernelInput( xi, xj, eps, geomKind )
	real(kreal) :: integralKernelInput
	real(kreal), intent(in) :: xi(3), xj(3), eps
	integer(kint), intent(in) :: geomKind	
	if ( geomKind == PLANAR_GEOM ) then
		integralKernelInput = ChordDistance(xi,xj)/eps
	elseif ( geomKind == SPHERE_GEOM ) then
		integralKernelInput = SphereDistance(xi,xj)/eps
	else	
		integralKernelInput = 0.0_kreal
	endif
end function

pure function bivariateDeltaKernel8( radialDist ) 
	real(kreal) :: bivariateDeltaKernel8
	real(kreal), intent(in) :: radialDist
	bivariateDeltaKernel8 = (4.0_kreal - 6.0_kreal * radialDist**2 + 2.0_kreal * radialDist**4 - &
			 radialDist**6 / 6.0_kreal) * exp(-radialDist * radialDist) / PI
end function

pure function bivariateLaplacianKernel8( radialDist )
	real(kreal) :: bivariateLaplacianKernel8
	real(kreal), intent(in) :: radialDist
	bivariateLaplacianKernel8 = (40.0_kreal - 40.0_kreal * radialDist**2 + 10.0_kreal * radialDist**4 - &
								 2.0_kreal * radialDist**6/3.0_kreal) * exp( -radialDist*radialDist) / PI
end function

pure function bivariateFirstDerivativeKernel8( radialDist)
	real(kreal) :: bivariateFirstDerivativeKernel8
	real(kreal), intent(in) :: radialDist
	bivariateFirstDerivativeKernel8 = (-20.0_kreal + 20.0_kreal * radialDist**2 - 5.0_kreal * radialDist**4 + &
									   radialDist**6/3.0_kreal) * exp( - radialDist*radialDist ) / PI
end function

pure function bivariateSecondDerivativeKernel8( radialDist )
	real(kreal) :: bivariateSecondDerivativeKernel8
	real(kreal), intent(in) :: radialDist
	bivariateSecondDerivativeKernel8 = (24.0_kreal -16.0_kreal * radialDist**2 + 2.0_kreal * radialDist**4) * &
				exp(- radialDist**2) / PI
end function

pure function bivariateMixedSecondDerivativeKernel8( radialDist, x, y )
	real(kreal) :: bivariateMixedSecondDerivativeKernel8
	real(kreal), intent(in) :: radialDist, x, y
	bivariateMixedSecondDerivativeKernel8 = x * y * (160.0_kreal - 120.0_kreal * radialDist**2  + &
		24.0_kreal * radialDist**4 - 4.0_kreal * radialDist**6/3.0_kreal )* exp(-radialDist**2)/PI
end function

pure function trivariateDeltaKernel8( radialDist )
	real(kreal) :: trivariateDeltaKernel8
	real(kreal), intent(in) :: radialDist
	trivariateDeltaKernel8 = (6.5625_kreal - 7.875 * radialDist**2 + 2.25_kreal * radialDist**4 - &
		radialDist**6/6.0_kreal) * exp( - radialDist**2 ) / 5.568327996831708_kreal
end function

pure function trivariateFirstDerivativeKernel8( radialDist, eps )
	real(kreal) :: trivariateFirstDerivativeKernel8
	real(kreal), intent(in) :: radialDist
	real(kreal), intent(in) :: eps
	trivariateFirstDerivativeKernel8 = (-7.875_kreal/eps**2 + 9.0_kreal * radialDist**2/eps**4 - radialDist**4/eps**4 - &
										13.125_kreal - 15.75_kreal * radialDist**2 + 4.5_kreal * radialDist**4 - &
										radialDist**6/3.0_kreal ) * exp(-radialDist**2)
end function

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine


end module 