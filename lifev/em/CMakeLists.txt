INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(EM)
SET(${PACKAGE_NAME}_VERSION "1.0")

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORY(solver)
SET(HEADERS ${HEADERS} 
            ${solver_HEADERS} 
            ${solver_electrophysiology_HEADERS}  
            ${solver_electrophysiology_ionicModels_HEADERS}
            ${solver_mechanics_HEADERS} 
            ${solver_mechanics_materials_HEADERS}
            ${solver_mechanics_materials_passive_HEADERS}
            ${solver_mechanics_materials_activeStress_HEADERS}
            ${solver_mechanics_materials_activeStrain_HEADERS}
            ${solver_activation_HEADERS}
            ${solver_activation_xBridgeModels_HEADERS}
            ${solver_activation_activeStressModels_HEADERS}
            ${solver_activation_activeStrainModels_HEADERS}
            ${solver_circulation_HEADERS})

SET(SOURCES ${SOURCES} 
            ${solver_SOURCES} 
            ${solver_electrophysiology_SOURCES} 
            ${solver_electrophysiology_ionicModels_SOURCES}
            ${solver_mechanics_SOURCES}
            ${solver_mechanics_materials_SOURCES}
            ${solver_mechanics_materials_passive_SOURCES}
            ${solver_mechanics_materials_activeStress_SOURCES}
            ${solver_mechanics_materials_activeStrain_SOURCES}
            ${solver_activation_SOURCES}
            ${solver_activation_xBridgeModels_SOURCES}
            ${solver_activation_activeStressModels_SOURCES}
            ${solver_activation_activeStrainModels_SOURCES}
            ${solver_circulation_SOURCES})

ADD_SUBDIRECTORY(util)
SET(HEADERS ${HEADERS} ${util_HEADERS})
SET(SOURCES ${SOURCES} ${util_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  lifevem
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)

TRIBITS_ADD_TEST_DIRECTORIES(testsuite)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
