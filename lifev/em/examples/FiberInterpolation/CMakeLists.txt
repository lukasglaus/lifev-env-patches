
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  FiberInterpolation
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_electrophysiology_FiberInterpolation_data
  CREATE_SYMLINK
  SOURCE_FILES data ParamList1.xml ParamList2.xml SolverParamList_rbf3d.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_FiberInterpolation
  SOURCE_FILES idealized.mesh ellipsoid.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/electrophysiology/data/mesh/
)
