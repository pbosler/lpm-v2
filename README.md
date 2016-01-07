# lpm-v2
Lagrangian Particle Methods for PDEs

LPM : Lagrangian Particle Method
=========

Code associated with 

> P. Bosler, J. Kent, R. Krasny, and C. Jablonowski. "A Lagrangian particle method with remeshing for tracer transport on the sphere,"
> submitted to _Journal of Computational Physics_, 2015.


> P. Bosler,  L. Wang,  C. Jablonowski, and R. Krasny.
>	"A Lagrangian particle/panel method for the barotropic vorticity equation on a rotating sphere," _Fluid Dynamics Research_,  46 : 031406, 2014.

> P. Bosler, "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis, University of Michigan, 2013.

Design / Use
=========

The code is organized in an object-oriented style using modern Fortran. 
Hence, there is a class (Fortran derived data type) that corresponds to every task that the code completes.
Generally, each data type has its own module file.@n
A C++ implementation is in development.

Base objects 
------------
LPM provides a small library of fundamental modules that are used by the rest of the code.

* The @ref NumberKinds module defines constants used by the rest of the code, including variable kind definitions and physical constants.
* The @ref OutputWriter class defines methods for formatting output either for the console or ASCII files.  
* The @ref Logger class provides objects and methods for outputting to the console various types of messages from the rest of the LPM code in a parallel computing environment.  
* The @ref SphereGeom and @ref PlaneGeom modules define geometric formulas for distance and areas. 
* The @ref MPISetup module provides a type and methods for distribution other objects across MPI ranks. 
* The @ref STDIntVector module provides a data type that mimics the C++ `std::vector<int>` class to dynamically allocate and resize integer arrays of rank 1

Particles and Fields
--------------------
The primary computing objects used by LPM are the Lagrangian particles, defined in the particles module, and the fields, defined in the field module. 
* @ref Particles provide the spatial discretization of a domain and have both physical and Lagrangian coordinates.  
* A @ref Field may be either a scalar or vector field defined on a set of particles; fields are used to store variables (e.g., vorticity, temperature, etc.) on a set of particles.  


Mesh objects
------------
While the dynamics and PDEs are solved by LPM in a ''meshfree'' context, we often find it helpful and efficient to maintain a mesh to organize our particles.   
Topologically two-dimensional meshes (e.g., spheres and planes) are currently implemented using triangles and quadrilaterals.  
* The @ref Edges module defines directed edges that connect an origin particle to a destination particle, and have a left face and a right face.
* The @ref Faces module defines polyhedral faces that have an associated list of edges and vertices.  
* The @ref PolyMesh2d module collects particles, edges, and faces, into a single mesh object.   
Meshes are initialized from a ''seed'' that is recursively divided until a desired spatial resolution is achieved.  

* @ref PolyMesh2d meshes may be adaptively refined using the @ref Refinement module.  
* Remesh/remap subroutines depend on a separate scattered data algorithm.  Current implementations include @ref BIVARRemesh for planar problems
and @ref SSRFPACKRemesh for spherical problems.

External libraries
-------------------
LPM code requires scattered data interpolation or approximation algorithms to perform its remesh/remap step.   
For planar applications, the bivariate quintic Hermite polynomials of the BIVAR pacakge (H. Akima, _ACM TOMS_, 1978) are used (bivar.f90).
For spherical applications, the cubic Hermite polynomials with exponential tension factors of the STRIPACK/SSRFPACK libraries (Renka, _ACM TOMS_, 1997) are used (stripack.f and ssrfpack.f).  

Other methods for scattered data approximation are under development.  

PDE types and Solver objects
--------------
PDEs are solved using a ''method of lines'' discretization. 
A mesh is connected to its relevant variables by a PDE type, e.g., an SWEMesh or a @ref SphereBVE BVEMesh object.
These objects combine a PolyMesh2d for the spatial discretization with the various Fields (e.g., vorticity, temperature) 
appropriate to the chosen PDE.

Following the spatial discretization given by the particles, the time-ODEs are integrated explicity with 4th order Runge-Kutta by a Solver
(e.g., @ref SphereBVESolver) object appropriate to a particular application, e.g., SWEPlaneSolver or SphereBVESolver.  

* Planar, inviscid, incompressible flow: See @ref PlanarIncompressible, @ref PlanarIncompressibleSolver, and CollidingDipoles.f90
* Rotating planar, inviscid, incompressible flow: See @ref BetaPlane, @ref BetaPlaneSolver
* Spherical Barotropic Vorticity Equation (BVE): See @ref SphereBVE, @ref SphereBVESolver, and BVESingleGaussianVortex.f90
* Rotating planar Shallow Water Equations (SWE): See @ref PlanarSWE, @ref SWEPlaneSolver, NitscheStrickland.f90, etc.

Application Drivers
--------
Finally, all of the above objects are assembled by a driver program that solves a specific problem.  
Within the driver program all initial conditions must be defined.

Run-time user input to driver programs is handled via namelist file, that is typically supplied as the first argument on the command line.
See below for more discussion.

It is within these driver programs that the timestepping loop and remeshing procedures are called.
See PlaneSWEGravityWaves.f90 or CollidingDipoles.f90 for examples.


Build / Install
================
This latest version of LPM is written with support for KitWare's [CMake](http://www.cmake.org) cross-platform makefile generator, to aid
portability and ease installation for new users.

LPM requires an [OpenMPI](http://www.open-mpi.org) distribution for parallel computing.

The environmental variables `CC`, `CXX`, and `FC` specify the compilers.  
For example, using the bash shell, these variables are set with

    export CC=mpicc
    export CXX=mpicxx
    export FC=mpifort
    
LPM is configured in the same shell by navigating to the lpm-v2 root directory and typing 
    mkdir build
    cd build
    cmake ..
This configures the build for the your compiler version and MPI.  Build the LPM software next by typing

    make
    
then

    make test
    
to run the unit tests.  

Running applications
====================

LPM uses Fortran's handy namelist utility to specify run-time input to programs.

To run an application driver, LPM executables may be run with the command

    mpirun -np 6 <LPMfile>.exe  <namelist file>

This would execute the specified .exe file using 6 MPI processes and the input variables defined in the specified namelist file.

Driver programs
----------------

Each application is defined in its own driver program.  
Examples include `CollidingDipoles.f90`, `BVESingleGaussianVortex.f90`, and `PlaneSWEGravityWaves.f90`. @n

Driver programs are organized into three phases:
1. __Initialize:__ The computing environment is defined and initialized using appropriate calls to MPI subroutines. @n
	User input is read from the appropriate namelist file. @n
	The initial set of particles and their initial mesh are defined. @n
	Initial conditions for each variable are defined on the mesh. @n
	Any initial adaptive refinement is completed. @n
	The initial conditions are output to data files. @n
2. __Run:__ The time step loop resides in the section of the driver. @n
	Within the timestep loop, an `if` statement determinces whether or not a particle set needs to be remeshed and/or remapped. @n
	Remeshing (if necessary) is completed. @n
	Time is advanced by one time step for each iteration of the loop. @n
	Time-dependent output is written to a data file. @n
3. __Finalize:__ Final output is written to data files. @n
	Performance and timing data are output to console. @n
	All LPM objects are deleted and all memory used by the driver is freed.

Within the `contains` section of a driver program, users may insert functions relevant to their application.  @n
In most cases, these functions must conform to the interfaces defined in @ref NumberKinds and @ref Refinement . @n
Additionally, since different applications will require different input, each driver requires its own version of the `ReadNamelistFile` subroutine.
