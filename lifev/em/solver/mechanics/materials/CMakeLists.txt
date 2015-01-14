ADD_SUBDIRECTORY(functions)
ADD_SUBDIRECTORY(passive_materials)
ADD_SUBDIRECTORY(active_strain_materials)
ADD_SUBDIRECTORY(active_stress_materials)

SET( solver_mechanics_materials_HEADERS
  solver/mechanics/materials/EMMaterialType.hpp
CACHE INTERNAL "")

SET( solver_mechanics_materials_SOURCES
# solver/mechanics/materials/EMMaterial.cpp
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_mechanics_materials_HEADERS})
  STRING(REPLACE "solver/mechanics/materials/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/mechanics/materials/"
  COMPONENT ${PACKAGE_NAME}
)
