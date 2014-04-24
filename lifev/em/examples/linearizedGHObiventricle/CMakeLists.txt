
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  linearizedGHObiventricle
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_linearizedGHObiventricle
   SOURCE_FILES ParamList.xml data
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(cartesian_cube8.mesh
#   SOURCE_FILES cartesian_cube8.mesh
#   SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/darcy/data/mesh/inria
# )

