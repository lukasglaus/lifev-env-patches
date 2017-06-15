
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
#basicExample
# example_InflationExtension
#  example_applyForcingTermOnDeformedConfiguration/
# example_CookMembrane
# example_creatingDamagedZone
# example_tractionWithSymmetry
#  example_principalTensionsInflationExtensions
# example_partitionMesh
# example_computationNorm
)

IF(LIVEF_HAS_HDF5)
ADD_SUBDIRECTORIES(
# example_applyForcingTermOnDeformedConfiguration/
)
ENDIF()
