
INCLUDE(TribitsAddExecutable)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
#  basicExample
#  basicExampleHO
#  basicExampleGHO
#  exampleGHO_withFiberRotation
#  exampleHO_withFiberRotation
#  exampleStatic_GAHO
#  example_ActiveStrain
#  exampleHO_ellipsoid
#  example_CalciumDependentActivation
#  example_laplacian
#  linearizedGHO
#  linearizedGHObiventricle
#  linearizedGHOrotation
#  linearizedGHOventricle
#  prestress
   example_EMMaterial
   example_HeartOharaRudy
#   example_passiveEMMaterial
   example_EMCoupling
   example_XBridge
   example_EMBiVent
   example_EMSiVent
)
