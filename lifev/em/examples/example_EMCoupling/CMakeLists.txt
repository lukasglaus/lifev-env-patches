
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(/usr/local/Cellar/eigen/3.2.4/include/eigen3)

TRIBITS_ADD_EXECUTABLE(
  EMCoupling
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_example_EMCoupling
   SOURCE_FILES ParamList.xml data dataP1 FiberDirectionP1.h5 FiberDirectionP2.h5
   SheetsDirectionP1.h5 SheetsDirectionP2.h5
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(cube4.mesh
#   SOURCE_FILES cube4.mesh
#   SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/data/mesh
# )

