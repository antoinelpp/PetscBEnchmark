cmake_minimum_required(VERSION 2.5)
enable_language (Fortran)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#-------------------------------------------------------------------
#Génération de la liste des fichiers sources
#-------------------------------------------------------------------
file(GLOB_RECURSE COMMON_SOURCES src psolver.c hsolver.c utilisis_mod.f90)
file(GLOB_RECURSE TARGET_SOURCES src main.f90)

message("Source1: ${COMMON_SOURCES}")
message("Source 2: ${TARGET_SOURCES}")


#------------------------------------------------------------------
#Find MPI compilers, with Intel compilers if found
#-------------------------------------------------------------------
#get the env compilator variable : "mpifort" or "mpiifort" (defined in .bashrc)
message("")
message("${Red}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
message("testing the host name${ColourReset}")

SITE_NAME(name)
message("You are compiling on:${name}")
if (${name} STREQUAL "lx-tavant.lab-lpp.local" OR
    ${name} STREQUAL "zoidberg.lab-lpp.local")
  message("${Green}You should use GNU fortran${ColourReset}")
elseif(  ${name} STREQUAL "nemo1" OR ${name} STREQUAL "nemo2")
  message("You should use Intel fortran")

else()
  message("the host is not know.")
endif()

message("")
message("${Yellow}~~~~~~~~ Detecting Compiler  ~~~~~~~~~~${ColourReset}")

set(MPI_Fortran_COMPILER $ENV{COMP})

message("Compilator wanted : ${MPI_Fortran_COMPILER}" )

if(${MPI_Fortran_COMPILER} STREQUAL "mpiifort")
  find_package(MPI REQUIRED)
  set( CMAKE_Fortran_COMPILER_ID  Intel)

else() # lest try MPIF90

  #By setting MPI_C_COMPILER before to the alias set by the "usr module", FindMPI will look for this specific compiler
  set( MPI_Fortran_COMPILER mpifort)
  set( MPI_C_COMPILER  mpicc)
  set( MPI_Fortran_COMPILER  mpifort)

  #find_package(MPI REQUIRED)

  set( CMAKE_Fortran_COMPILER_ID  GNU)

endif()

#Set CMAKE variables
set( CMAKE_C_COMPILER  ${MPI_C_COMPILER})
set( CMAKE_Fortran_COMPILER  ${MPI_Fortran_COMPILER})

#For small debug : MPIversion will be either mpifort or mpiifort
get_filename_component(MPIversion ${MPI_Fortran_COMPILER} NAME)

message("FORTRAN COMPILER =     ${MPI_Fortran_COMPILER}")

#include and link
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_COMPILE_FLAGS ${CMAKE_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})

# Choix des options de compilations
#------------------------------------------------------------------

# DEBUG: if you set the environement variable DEBUG to True, the code is compile with debug flags.
# Else, we optimise by default.
# example : `$ DEBUG=True cmake ..` will enable the debuger
set(DEBUG $ENV{DEBUG})
set(DEBUG true)

if(DEBUG)
  message("You choose to ${BoldRed}debug${ColourReset}")
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set( CMAKE_Fortran_FLAGS "-check all -CB -g -traceback" )
    set( CMAKE_C_FLAGS "-g -rdynamic")
  else()
    set( CMAKE_Fortran_FLAGS "-g -fbacktrace" )
    set( CMAKE_C_FLAGS "-g -rdynamic")

  endif()

else(DEBUG)
  message("You choose to ${BoldRed}Optimize${ColourReset} !")
    set( CMAKE_Fortran_FLAGS "-O2" )
    set( CMAKE_C_FLAGS "-O2")

endif(DEBUG)
message("Compilation flags : ${CMAKE_Fortran_FLAGS}")

#------------------------------------------------------------------
#SET HDF5
#------------------------------------------------------------------
message("--------------------   HDF5     --------------------------")
set(HDF5_DIR $ENV{HDF5_DIR})
set(HDF5_INCLUDE_DIRS $ENV{HDF5_INCLUDE})
if(EXISTS ${HDF5_DIR})
  message("HDF5_DIR = ${HDF5_DIR}")
  set(HDF5_LIB_DIR  ${HDF5_DIR}/lib)
  if(NOT EXISTS ${HDF5_INCLUDE_DIRS})
      set(HDF5_INCLUDE_DIRS ${HDF5_DIR}/include)
  endif()
else()
  message("${Red}HDF5_DIR in not SET. We will use HDF5_INCLUDE and HDF5_LIB${ColourReset}")
  set(HDF5_INCLUDE_DIRS $ENV{HDF5_INCLUDE})
  set(HDF5_LIB_DIR  $ENV{HDF5_LIB})
endif()

set(HDF5_Fortran_LIBRARIES ${HDF5_LIB_DIR}/libhdf5.so ${HDF5_LIB_DIR}/libhdf5_fortran.so  )

MESSAGE("HDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
MESSAGE("HDF5_LIB_DIR=${HDF5_LIB_DIR}")
MESSAGE("HDF5_Fortran_LIBRARIES=${HDF5_Fortran_LIBRARIES}")

#------------------------------------------------------------------
#FInd HYPRE
#------------------------------------------------------------------
message("--------------------   HYPRE     --------------------------")
set(HYPRE_DIR $ENV{HYPRE_DIR})
set(HYPRE_INCLUDE_DIRS $ENV{HYPRE_INCLUDE})
if(EXISTS ${HYPRE_DIR})
  message("HYPRE_DIR = ${HYPRE_DIR}")
  set(HYPRE_LIB_DIR  ${HYPRE_DIR}/lib)
  if(NOT EXISTS ${HYPRE_INCLUDE_DIRS})
    set(HYPRE_INCLUDE_DIRS ${HYPRE_DIR}/include)
  endif()
else()
  message("${Red}HYPRE_DIR in not SET. We will use HYPRE_INCLUDE and HYPRE_LIB${ColourReset}")
  set(HYPRE_INCLUDE_DIRS $ENV{HYPRE_INCLUDE})
  set(HYPRE_LIB_DIR  $ENV{HYPRE_LIB})
endif()

MESSAGE("HYPRE_INCLUDE_DIRS=${HYPRE_INCLUDE_DIRS}")
MESSAGE("HYPRE_LIB_DIR=${HYPRE_LIB_DIR}")
#------------------------------------------------------------------
#FInd PetSc
#------------------------------------------------------------------
message("--------------------   HYPRE     --------------------------")
set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_INCLUDE_DIRS $ENV{PETSC_INCLUDE})

if(EXISTS ${PETSC_DIR})
  message("PETSC_DIR = ${PETSC_DIR}")
  set(PETSC_LIB_DIR  ${PETSC_DIR}/lib)
  if(NOT EXISTS ${PETSC_INCLUDE_DIRS})
  	set(PETSC_INCLUDE_DIRS ${PETSC_DIR}/include)
  endif()
else()
  message("${Red}PETSC_DIR in not SET. We will use PETSC_INCLUDE and PETSC_LIB${ColourReset}")
  set(PETSC_INCLUDE_DIRS $ENV{PETSC_INCLUDE})
  set(PETSC_LIB_DIR  $ENV{PETSC_LIB})
endif()

set(PETSC_LIBRARY ${PETSC_LIB_DIR}/libpetsc.so )

MESSAGE("PETSC_INCLUDE_DIRS = ${PETSC_INCLUDE_DIRS}")
MESSAGE("PETSC_LIB_DIR = ${PETSC_LIB_DIR}")
message("PETSC Lib is = ${PETSC_LIBRARY}")

#------------------------------------------------------------------
#------------------------------------------------------------------
message("rep source du projet: ${LPPic2D_SOURCE_DIR}")
message("rep compilation du projet: ${LPPic2D_BINARY_DIR}")
message("")
message("--------------------   Links     --------------------------")
message("")
#Inclusion des biblios
include_directories(
              ${MPI_INCLUDE_PATH}
              ${HDF5_INCLUDE_DIRS}
              ${HYPRE_INCLUDE_DIRS}
              ${PETSC_INCLUDE_DIRS}
)

link_directories(
              ${PETSC_LIB_DIR}
              ${HDF5_LIB_DIR}
              ${HYPRE_LIB_DIR}
)

#Déclaration de l'exécutable
add_executable(lppic ${COMMON_SOURCES} ${TARGET_SOURCES})

#------------------------------------------------------------------
#------------------------------------------------------------------
#Configuration de l'édition de liens
target_compile_options(lppic PRIVATE -DHAVE_CONFIG_H -DHYPRE_TIMING)
target_link_libraries(

        lppic

        ${CMAKE_Fortran_FLAGS}

        ${CMAKE_C_FLAGS}

        ${MPI_LIBRARIES}

        ${PETSC_LIBRARY}

        ${HDF5_Fortran_LIBRARIES}

        HYPRE
)
