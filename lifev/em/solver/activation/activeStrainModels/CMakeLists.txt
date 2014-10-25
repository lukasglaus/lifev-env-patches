SET( solver_activation_activeStressModels_HEADERS
  solver/activation/activeStrainModels/ActiveStrainRossiModel14.hpp
  solver/activation/activeStrainModels/ActiveStrainActivation.hpp
CACHE INTERNAL "")

SET( solver_activation_activeStressModels_SOURCES
  solver/activation/activeStrainModels/ActiveStrainRossiModel14.cpp
 CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_acitvation_activeStrainModels_HEADERS})
  STRING(REPLACE "solver/activation/activeStrainModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/activation/activeStrainModels/"
  COMPONENT ${PACKAGE_NAME}
)
