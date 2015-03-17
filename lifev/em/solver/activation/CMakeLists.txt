ADD_SUBDIRECTORIES(
activeStressModels
activeStrainModels
xBridgeModels
)

SET( solver_activation_HEADERS
  solver/activation/Activation.hpp
  solver/activation/ActivationModelsList.hpp
CACHE INTERNAL "")

SET( solver_activation_SOURCES
# solver/mechanics/materials/EMMaterial.cpp
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_activation_HEADERS})
  STRING(REPLACE "solver/activation/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/activation/"
  COMPONENT ${PACKAGE_NAME}
)
