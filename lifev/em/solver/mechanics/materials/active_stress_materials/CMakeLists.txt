SET( solver_mechanics_materials_activeStress_HEADERS
 solver/mechanics/materials/active_stress_materials/SimpleFibersActiveStress.hpp
 CACHE INTERNAL "")

SET( solver_mechanics_materials_activeStress_SOURCES
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_mechanics_materials_activeStress_HEADERS})
  STRING(REPLACE "solver/mechanics/materials/active_stress_materials/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/mechanics/materials/active_stress_materials/"
  COMPONENT ${PACKAGE_NAME}
)
