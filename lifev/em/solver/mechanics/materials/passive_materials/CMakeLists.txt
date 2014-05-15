SET( solver_mechanics_materials_passive_HEADERS
 solver/mechanics/materials/passive_materials/PassiveIsotropicExponential.hpp
 solver/mechanics/materials/passive_materials/PassiveLinearizedNeoHookean.hpp
 solver/mechanics/materials/passive_materials/PassiveNeoHookean.hpp
 CACHE INTERNAL "")

SET( solver_mechanics_materials_passive_SOURCES
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_mechanics_materials_passive_HEADERS})
  STRING(REPLACE "solver/mechanics/materials/passive_materials/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/mechanics/materials/passive_materials/"
  COMPONENT ${PACKAGE_NAME}
)
