INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
#  EMBenchmark
#  EMCell
#  EMVentricle
#  EMMyocyte
	test_passiveEMMaterial
	test_activeNeoHook
	test_OneWayCoupling
	test_fixedMonodomain
	test_benchmarkBar	
	test_benchmarkStaticActiveStress
	test_benchmarkIsotropicVentricle
	test_HDF5toVTK
	test_EMSolver
)
