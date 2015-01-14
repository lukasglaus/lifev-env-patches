ADD_SUBDIRECTORY(IonicModels)

SET( solver_electrophysiology_HEADERS
	solver/electrophysiology/EMMonodomainSolver.hpp
    solver/electrophysiology/IonicModelsList.hpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_electrophysiology_HEADERS})
  STRING(REPLACE "solver/electrophysiology/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/electrophysiology"
  COMPONENT ${PACKAGE_NAME}
)
