
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE(
  test_generateFibers
  SOURCES main_fibers.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  test_StaticActiveStress
  SOURCES main_benchmark.cpp
  )


 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_staticActiveStress
   SOURCE_FILES ParamList.xml data_fibers data_benchmark
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ellipsoid_staticActiveStress
  SOURCE_FILES ellipsoid_5mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/data/mesh/
)

