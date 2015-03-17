SET( solver_activation_xBridgeModels_HEADERS
  solver/activation/xBridgeModels/xBridge4SM.hpp
CACHE INTERNAL "")

SET( solver_activation_xBridgeModels_SOURCES
  solver/activation/xBridgeModels/ActiveStrainRossiModel14.cpp
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_acitvation_xBridgeModels_HEADERS})
  STRING(REPLACE "solver/activation/xBridgeModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/activation/xBridgeModels/"
  COMPONENT ${PACKAGE_NAME}
)
