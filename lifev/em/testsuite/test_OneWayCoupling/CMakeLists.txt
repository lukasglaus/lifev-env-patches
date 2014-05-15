
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE(
  test_OneWayCoupling
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_OneWayCoupling
   SOURCE_FILES ParamList.xml data
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

TRIBITS_COPY_FILES_TO_BINARY_DIR(cube4_test_OneWayCoupling
  SOURCE_FILES cube4.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/data/mesh/
)

