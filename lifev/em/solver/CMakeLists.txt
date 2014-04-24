ADD_SUBDIRECTORIES(
	electrophysiology
	Mechanics
)

SET(solver_HEADERS
#  solver/EMActiveStrainSolver.hpp
#  solver/EMEvaluate.hpp
  solver/EMSolver.hpp
#  solver/EMETAFunctors.hpp
#  solver/EMGeneralizedActiveHolzapfelOgdenMaterial.hpp
#  solver/EMActiveStructuralConstitutiveLaw.hpp
#  solver/EMNeoHookeanActivatedMaterial.hpp
#  solver/EMHolzapfelOgdenMaterial.hpp
CACHE INTERNAL "")

SET(solver_SOURCES
  solver/EMSolver.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_HEADERS})
  STRING(REPLACE "solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/em/solver"
  COMPONENT ${PACKAGE_NAME}
)
