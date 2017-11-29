SET( solver_mechanics_materials_activeStrain_HEADERS
 solver/mechanics/materials/active_strain_materials/ActiveNeoHookean.hpp
solver/mechanics/materials/active_strain_materials/ActiveIsotropicExponential.hpp
solver/mechanics/materials/active_strain_materials/ActiveHolzapfelOgden.hpp
solver/mechanics/materials/active_strain_materials/ActiveHolzapfelOgdenNIIsotropicSplitting.hpp
 CACHE INTERNAL "")

SET( solver_mechanics_materials_activeStrain_SOURCES
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_mechanics_materials_activeStrain_HEADERS})
  STRING(REPLACE "solver/mechanics/materials/active_strain_materials/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/mechanics/materials/active_strain_materials/"
  COMPONENT ${PACKAGE_NAME}
)
