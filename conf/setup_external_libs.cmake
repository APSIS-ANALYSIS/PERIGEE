# This File Setup the value of EXTRA_LINK_LIBS and the INCLUDE_DIRS

MESSAGE(STATUS "=================== External library setup ===================")
# ===============================================
# SET VTK
# ===============================================
IF(DEFINED VTK_DIR)
  MESSAGE(STATUS "VTK_DIR = ${VTK_DIR}")
ELSE(DEFINED VTK_DIR)
  MESSAGE(FATAL_ERROR "VTK_DIR is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED VTK_DIR)

If(DEFINED VTK_VERSION)
  MESSAGE(STATUS "VTK_VERSION = ${VTK_VERSION}")
ELSE(DEFINED VTK_VERSION)
  MESSAGE(FATAL_ERROR "VTK_VERSION is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED VTK_VERSION)

IF(DEFINED VTK_link_lib)
  MESSAGE(STATUS "VTK_link_lib=${VTK_link_lib}")
ELSE(DEFINED VTK_link_lib)
  MESSAGE(FATAL_ERROR "VTK_link_lib is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED VTK_link_lib)

SET(VTK_INCLUDE_DIRS ${VTK_INCLUDE_DIRS}
  ${VTK_DIR}/include/${VTK_VERSION})


IF(${VTK_VERSION} MATCHES "vtk-5.10")
  SET(VTK_LINK_DIRS ${VTK_LINK_DIRS} ${VTK_DIR}/lib/${VTK_VERSION})
ELSE(${VTK_VERSION} MATCHES "vtk-5.10")
  SET(VTK_LINK_DIRS ${VTK_LINK_DIRS} ${VTK_DIR}/lib)
ENDIF(${VTK_VERSION} MATCHES "vtk-5.10")

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_link_lib} )

# ===============================================
# FINISH VTK SETUP
# ===============================================

# ===============================================
# SET PETSc library 
# ===============================================
IF(DEFINED PETSC_DIR)
  MESSAGE(STATUS "PETSc directory = ${PETSC_DIR}")
ELSE(DEFINED PETSC_DIR)
  MESSAGE(FATAL_ERROR "PETSC_DIR is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_DIR)

IF(DEFINED PETSC_ARCH)
  MESSAGE(STATUS "PETSc arch = ${PETSC_ARCH}")
ELSE(DEFINED PETSC_ARCH)
  MESSAGE(FATAL_ERROR "PETSC_ARCH is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_ARCH)

IF(DEFINED PETSC_LIBRARY_DIRS)
  MESSAGE(STATUS "PETSC_LIBRARY_DIRS = ${PETSC_LIBRARY_DIRS}")
ELSE(DEFINED PETSC_LIBRARY_DIRS)
  MESSAGE(FATAL_ERROR "PETSC_LIBRARY_DIRS is NOT defined! Check ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_LIBRARY_DIRS)

IF(DEFINED PETSC_link_lib)
  MESSAGE(STATUS "PETSC_link_lib = ${PETSC_link_lib}")
ELSE(DEFINED PETSC_link_lib)
  MESSAGE(FATAL_ERROR "PETSC_link_lib is NOT defined! Check
  ~/PETSc_VTK.cmake")
ENDIF(DEFINED PETSC_link_lib)

# Add the PETSc libs to the lib-link-list
SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_link_lib})

SET(PETSC_INCLUDE_DIRS ${PETSC_DIR}/include
  ${PETSC_DIR}/${PETSC_ARCH}/include)
# ===============================================
# Finish PETSc Library Setup
# ===============================================

# ===============================================
# SET METIS library
# ===============================================
IF(DEFINED METIS_DIR)
  MESSAGE(STATUS "METIS DIR = ${METIS_DIR}")
ELSE(DEFINED METIS_DIR)
  MESSAGE(FATAL_ERROR "METIS_DIR is NOT defined! check ~/PETSc_VTK.cmake")
ENDIF(DEFINED METIS_DIR)

SET(METIS_LIBRARY_DIRS ${METIS_DIR}/lib)
SET(METIS_INCLUDE_DIRS ${METIS_DIR}/include)
SET(METIS_LINK_LIBS metis)

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LINK_LIBS})
# ===============================================
# FINISH METIS setup
# ===============================================

# ===============================================
# SET HDF5 library
# ===============================================
IF(DEFINED HDF5_DIR)
  MESSAGE(STATUS "HDF5 DIR = ${HDF5_DIR}")
ELSE(DEFINED HDF5_DIR)
  MESSAGE(FATAL_ERROR "HDF5_DIR is NOT defined! check ~/PETSc_VTK.cmake")
ENDIF(DEFINED HDF5_DIR)

SET(HDF5_LIBRARY_DIRS ${HDF5_DIR}/lib)
SET(HDF5_INCLUDE_DIRS ${HDF5_DIR}/include)
SET(HDF5_LINK_LIBS hdf5)

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LINK_LIBS})
# ===============================================
# FINISH HDF5 setup
# ===============================================

# ===============================================
# Finish Slepc setup
# ===============================================
# INCLUDE LIBRARY HEADER FILES
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${PETSC_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${METIS_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )

# LINK TO LIBRARY LIBS
LINK_DIRECTORIES(${VTK_LINK_DIRS})
LINK_DIRECTORIES(${PETSC_LIBRARY_DIRS})
LINK_DIRECTORIES(${METIS_LIBRARY_DIRS})
LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})

# EOF
