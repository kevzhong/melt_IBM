# melt_IBM
This is a repository for a DNS code which simulates the Lagrangian melting of a solid object. The journal paper detailing this work is found at:

* Zhong, K., Howland, C. J., Lohse, D., and Verzicco, R. (2025). A front-tracking immersed-boundary framework for simulating Lagrangian melting problems. J. Comput. Phys., page 113762.

There is also a WIP documentation for this repository found [here](https://kevzhong.github.io/melt_IBM/).

# Dependencies
- MPI Fortran compiler
- HDF5
- FFTW3
- LAPACK, BLAS, etc

Fortran compiler can be GNU (`gfortran`) or Intel, but it is important that the HDF5 library was compiled using the specified Fortran compiler.

For MacOS, dependendices can be installed simply:

`brew install fftw hdf5-mpi open-mpi openblas lapack gcc`

# Usage
Simulation parameters are specified in `bou.in` and `part.in`. HIT parameters specified in `HITForcing.in`. 

# Key features
- Incompressible Navier--Stokes + temperature equation, with optional coupling through thermal Boussinesq driving
- Second-order staggered, finite difference + low-storage RK3 substepping
- Moving-least-squares immersed boundary method treatment: immersed solid is a triangulated surface, which is read-in as a `.gts` file
- Fluid-solid interaction governed by Newton--Euler equations
- Melting governed by the Stefan condition
- Triangulated solid geometry is dynamically remeshed as it melts
- MPI parallelisation is the slab decomposition

# Main restrictions
- Single solid object (no collision treatment)
- Periodic domain
- Matched thermal diffusivities between solid and liquid phases
- Matched densities between the phases