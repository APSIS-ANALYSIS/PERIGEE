# This File Setup the value of EXTRA_LINK_LIBS and the INCLUDE_DIRS

# ===============================================
# SET VTK
# ===============================================
SET(VTK_INCLUDE_DIRS ${VTK_INCLUDE_DIRS}
  ${VTK_DIR}/include/${VTK_VERSION})

IF(${VTK_VERSION} MATCHES "vtk-5.10")
  SET(VTK_LINK_DIRS ${VTK_LINK_DIRS} ${VTK_DIR}/lib/${VTK_VERSION})
ELSE(${VTK_VERSION} MATCHES "vtk-5.10")
  SET(VTK_LINK_DIRS ${VTK_LINK_DIRS} ${VTK_DIR}/lib)
ENDIF(${VTK_VERSION} MATCHES "vtk-5.10")

SET(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_link_lib} )

# ===============================================
# SET HDF5 library
# ===============================================
# INCLUDE LIBRARY HEADER FILES
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

# LINK TO LIBRARY LIBS
LINK_DIRECTORIES(${VTK_LINK_DIRS})

# EOF
