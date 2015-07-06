
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  EMCoupling
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_example_EMCoupling
   SOURCE_FILES ParamList.xml dataP1BiVent dataP2BiVent dataP1SiVent inputfile singleVentricle.mesh FiberDirection.h5  SheetsDirection.h5
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )