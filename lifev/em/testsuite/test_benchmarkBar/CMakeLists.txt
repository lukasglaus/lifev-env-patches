
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE(
  test_benchmarkBar
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(datatest_benchmarkBar
   SOURCE_FILES ParamList.xml data
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
   CREATE_SYMLINK
 )

TRIBITS_COPY_FILES_TO_BINARY_DIR(cube4_benchmarkBar
  SOURCE_FILES bar_05mm.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/data/mesh/
)

