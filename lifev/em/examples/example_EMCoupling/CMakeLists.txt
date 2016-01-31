
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  EMCoupling
  SOURCES main.cpp
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(data_example_EMCoupling
   SOURCE_FILES ParamList.xml data dataP1BiVent dataP2BiVent dataP1SiVent dataP2SiVent inputfile singleVentricle.mesh FiberDirection.h5  SheetsDirection.h5 siVentP1/FiberDirection.h5 siVentP1/SheetsDirection.h5 siVentP2/FiberDirection.h5 siVentP2/SheetsDirection.h5 xVentSingle.mesh singleVentConv/FiberDirection.h5 singleVentConv/SheetsDirection.h5 singleVentConv/singleVentConv.mesh
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
 )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(biVentCoarse
#  CREATE_SYMLINK
#  SOURCE_FILES biVentCoarse
#  SOURCE_DIR ${CMAKE_SOURCE_DIR}/../lifev-env/lifev-em-build/lifev/electrophysiology/testsuite/test_fibersHeart/
#)
