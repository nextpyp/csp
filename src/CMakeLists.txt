PROJECT( csp )

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

### START: ADDED BY ABARTE ###
#SET(ITK_DIR "${HOME}/data/ETTK-workspace/external/InsightToolkit-3.20.0")
#SET(VTK_DIR "${HOME}//data/ETTK-workspace/external/InsightToolkit-3.20.0")
#SET(ITK_DIR "/hpc/group/bartesaghilab/jwm14/ETTK/external/InsightToolkit-4.2.1/build")
#SET(VTK_DIR "/hpc/group/bartesaghilab/jwm14/ETTK/external/VTK5.10.1")
	
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -m32 -DBZDEBUG")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -no-multibyte-chars -fabi-version=11 -wd10237 -O3 -no-prec-div -no-prec-sqrt -pthread -static-intel -DMKL -mkl=sequential -g")

# Use mpi c++ compiler
# SET(CMAKE_CXX_COMPILER /opt/apps/rhel7/mpich-3.2.1/bin/mpic++)
#SET(MPI_DIR /usr/local/mpich2_i386)
# SET(MPI_DIR /opt/apps/rhel7/mpich-3.2.1)
# SET(MPI_LIBRARIES mpichcxx mpich pthread rt)

#SET(MPI_INCLUDE_DIR ${MPI_DIR}/include)
#SET(MPI_LIBRARY_DIR ${MPI_DIR}/lib)
# Configure Blitz
# FIND_PATH(BLITZ_DIR blitz/array.h ${CMAKE_CURRENT_SOURCE_DIR}/external/blitz-0.9)
# SET(BLITZ_DIR /hpc/group/bartesaghilab/jwm14/ETTK/external/blitz-0.9)
SET(BLITZ_INCLUDE_DIR ${BLITZ_DIR})
SET(BLITZ_LIBRARY_DIR ${BLITZ_DIR}/lib)

### END:   ADDED BY ABARTE ###

# Locate ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# Locate VTK.
FIND_PACKAGE ( VTK )
IF ( VTK_FOUND )
     INCLUDE(${VTK_USE_FILE})
ENDIF( VTK_FOUND )

# FIND_PACKAGE(OpenMP)
# if (OPENMP_FOUND)
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# else(OPENMP_FOUND)
#   MESSAGE( FATAL_ERROR "Compiler must support OpenMP")  
# endif(OPENMP_FOUND)

SET (
  csp
  flrSPAEulerAngles.cxx
  flrClasses.cxx
  flrCostFunctions.cxx
  flrReadParametersFile.cxx
  flrMain.cxx
)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
    # ${MPI_INCLUDE_DIR}			# ADDED BY ABARTE
  ${BLITZ_INCLUDE_DIR}	     	# ADDED BY ABARTE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# ADDED BY ABARTE
LINK_DIRECTORIES(
    # ${MPI_LIBRARY_DIR}
  ${BLITZ_LIBRARY_DIR}
  ${ITK_LIBRARY_DIRS}
)

# ADDED BY ABARTE
LINK_LIBRARIES(
    #${MPI_LIBRARIES}
    ${BLITZ_LIBRARIES}
    ${CISTEM_LIB}
)



ADD_EXECUTABLE( csp ${csp} ${UISrcs} )
ADD_EXECUTABLE( gaussian_distribution gaussian_distribution.cpp )
ADD_EXECUTABLE( SPA_euler_angles SPA_euler_angles.cpp flrSPAEulerAngles.cxx )
ADD_EXECUTABLE( SPA_statistics SPA_statistics.cpp )

SET (
  external_libraries
  vtkRendering
  vtkGraphics
  vtkIO
  vtkCommon
  vtkImaging
  vtkFiltering
  ${ITK_LIBRARIES}
)

# MESSAGE( "${external_libraries}" )

TARGET_LINK_LIBRARIES( csp ${external_libraries} )
# TARGET_LINK_LIBRARIES( SPA_euler_angles ${external_libraries} )

