INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  benchmark
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_EMbenchmark_data
  CREATE_SYMLINK
  SOURCE_FILES data ParamList.xml stimulus
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_EMbenchmark
  SOURCE_FILES benchmark_05mm.mesh benchmark_solid.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/data/mesh/inria/
)
