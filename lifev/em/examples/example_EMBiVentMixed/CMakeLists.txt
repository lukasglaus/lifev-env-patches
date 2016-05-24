
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  BiVentMixed
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_example_EMBiVentMixed
   SOURCE_FILES ParamList.xml data circulation
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

TRIBITS_COPY_FILES_TO_BINARY_DIR(biVentCoarsePatchesP2
CREATE_SYMLINK
SOURCE_FILES biVentCoarsePatchesP2
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/examples/example_EMBiVentMixed
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(biVentCoarsePatches
CREATE_SYMLINK
SOURCE_FILES biVentCoarsePatches
SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/examples/example_EMBiVentMixed
)
