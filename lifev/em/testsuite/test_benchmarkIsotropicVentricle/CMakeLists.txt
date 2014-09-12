
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE(
  test_benchmarkIsotropicVentricle
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(datatest_benchmarkIsotropicVentricle
   SOURCE_FILES ParamList.xml data
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ellipsoid_benchmarkIsotropicVentricle
  SOURCE_FILES ellipsoid_5mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/data/mesh/
)

