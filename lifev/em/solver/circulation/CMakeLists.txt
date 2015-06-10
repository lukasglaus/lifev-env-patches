INCLUDE_DIRECTORIES(/usr/local/Cellar/eigen/3.2.4/include/eigen3)

SET( solver_circulation_HEADERS
	solver/circulation/Circulation.hpp
CACHE INTERNAL "")

SET(LOCAL_HEADERS)
FOREACH(INC ${solver_circulation_HEADERS})
  STRING(REPLACE "solver/circulation/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver/circulation"
  COMPONENT ${PACKAGE_NAME}
)
