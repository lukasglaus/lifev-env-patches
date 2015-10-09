
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  BiVent
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_example_EMBiVent
   SOURCE_FILES ParamList.xml data inputfile
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

TRIBITS_COPY_FILES_TO_BINARY_DIR(biVentMedium
  CREATE_SYMLINK
  SOURCE_FILES biVentMedium
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/em/examples/example_EMBiVent/
)
