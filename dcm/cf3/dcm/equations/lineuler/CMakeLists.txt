# coolfluid_dcm_lineuler

list( APPEND coolfluid_dcm_equations_lineuler_files
  BCExtrapolation.hpp
  BCExtrapolation.cpp
  BCFarfield.hpp
  BCFarfield.cpp
  BCMirror.hpp
  BCMirror.cpp
  ComputeCharacteristicVariablesUniform2D.hpp
  ComputeCharacteristicVariablesUniform2D.cpp
  LibLinEuler.hpp
  LibLinEuler.cpp
  LinEuler2D.hpp
  LinEuler2D.cpp
  SourceMonopoleUniform.hpp
  SourceMonopoleUniform.cpp
  RightHandSide2D.hpp
  RightHandSide2D.cpp
  InitAcousticPulse.hpp
  InitAcousticPulse.cpp
  InitAcousticVorticityPulse.hpp
  InitAcousticVorticityPulse.cpp
  Terms.hpp
)

coolfluid_mark_not_orphan(
  BCWallNoSlip2D.hpp
  BCWallNoSlip2D.cpp
  BCWallNonUniformMeanflow2D.hpp
  BCWallNonUniformMeanflow2D.cpp
  ConvectionNonUniformMeanflow2D.hpp
  ConvectionNonUniformMeanflow2D.cpp
  Convection3D.hpp
  Convection3D.cpp
  FixMeanFlowBoundaryLayer.hpp
  FixMeanFlowBoundaryLayer.cpp
  SourceMeanflow2D.hpp
  SourceMeanflow2D.cpp
  SourceDipole2D.hpp
  SourceDipole2D.cpp
  SourceQuadrupole2D.hpp
  SourceQuadrupole2D.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_equations_lineuler 
                        SOURCES ${coolfluid_dcm_equations_lineuler_files}
                        LIBS    coolfluid_solver
                                coolfluid_dcm_core 
                                coolfluid_physics_lineuler
                        PLUGINS Physics )

