SET( solver_circulation_HEADERS
	solver/circulation/Circulation.hpp
    solver/circulation/CirculationAssembler.hpp
    solver/circulation/CirculationBCHandler.hpp
    solver/circulation/CirculationCoupling.hpp
    solver/circulation/CirculationDofHandler.hpp
    solver/circulation/CirculationElementFactory.hpp
    solver/circulation/CirculationGridView.hpp
    solver/circulation/CirculationIO.hpp
    solver/circulation/CirculationTimeIntegrator.hpp
    solver/circulation/CirculationVolumeIntegrator.hpp
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
