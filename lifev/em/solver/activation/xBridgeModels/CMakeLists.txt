SET( solver_activation_xBridgeModels_HEADERS
  solver/activation/XBridgeModels/XBridge4SM.hpp
CACHE INTERNAL "")

SET( solver_activation_xBridgeModels_SOURCES
  solver/activation/XBridgeModels/XBridge4SM.cpp
 CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_activation_xBridgeModels_HEADERS})
  STRING(REPLACE "solver/activation/XBridgeModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/activation/XBridgeModels/"
  COMPONENT ${PACKAGE_NAME}
)
