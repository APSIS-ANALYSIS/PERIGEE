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
#INCLUDE_DIRECTORIES( ${METIS_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )

# LINK TO LIBRARY LIBS
LINK_DIRECTORIES(${VTK_LINK_DIRS})
LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})

# EOF
