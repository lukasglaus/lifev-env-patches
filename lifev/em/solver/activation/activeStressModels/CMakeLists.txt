SET( solver_activation_activeStressModels_HEADERS
  solver/activation/activeStressModels/ActiveStressRossiModel14.hpp
  solver/activation/activeStressModels/ActiveStressNashPanfilovModel04.hpp
  solver/activation/activeStressModels/ActiveStressActivation.hpp
CACHE INTERNAL "")

SET( solver_activation_activeStressModels_SOURCES
  solver/activation/activeStressModels/ActiveStressRossiModel14.cpp
  solver/activation/activeStressModels/ActiveStressNashPanfilovModel04.cpp
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_activation_activeStressModels_HEADERS})
  STRING(REPLACE "solver/activation/activeStressModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/activation/activeStressModels/"
  COMPONENT ${PACKAGE_NAME}
)
