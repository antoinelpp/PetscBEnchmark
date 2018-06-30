# PetscBenchmarks

This Project propose a simple linear algebra problem to be solve with differents KSP and PC.

## Objective

It may be difficult to chose an appropriate Preconditionneur and solver for linear algrebra problems.

Indeed, it is really case dependant.

The obsective of this Benchmark is to help me chose a solver and preconditioneur easily, depending of the system size, conditions, and architecture (herdware).
Right know, I'm using only `PetsC`, but I am planning to add `Hypre` as well.

## The problem

The problem is a simple Poisson equation : $ \Delta \phi = \rho$

We are using Derichelt conditions in one direction, and Periodic conditions in the other.
The boundary conditions can be changed.

# Code architecture
This is a bit Tricky..
I am using the solver library from a Fortran code.
However, the Fortran API is not always presents depending of the hardware.
So, I use `*.c` files to call the solver, and I use a Fortran file for the main.

**But**, In order to ease the Benchmarking, I need to passe commandline arguments to the c files so that Petsc can access it !

## How to use it

First compile the code with the `cmake`. Be sure to define the environment variables like `Petsc_Dir` so that the library can be linked.

The execute the code with different command lines.
The file [`test\run.sh`](test\run.sh) gives an example.
