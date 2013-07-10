# coolfluid_dcm_euler

list( APPEND coolfluid_dcm_equations_euler_files
  LibEuler.hpp
  LibEuler.cpp
  BCExtrapolation.hpp
  BCExtrapolation.cpp
  BCFarfield.hpp
  BCFarfield.cpp
  BCMirror.hpp
  BCMirror.cpp
  BCPressureOutlet.hpp
  BCPressureOutlet.cpp
  BCPressureVelocityInlet.hpp
  BCPressureVelocityInlet.cpp
  BCSubsonicInlet.hpp
  BCSubsonicInlet.cpp
  Euler1D.hpp
  Euler1D.cpp
  Euler2D.hpp
  Euler2D.cpp
  RightHandSide1D.hpp
  RightHandSide1D.cpp
  RightHandSide2D.hpp
  RightHandSide2D.cpp
  Terms.hpp
  RiemannSolvers.hpp
  RiemannSolvers.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_equations_euler
                        SOURCES ${coolfluid_dcm_equations_euler_files}
                        LIBS    coolfluid_solver coolfluid_dcm_core coolfluid_physics_euler
                        PLUGINS Physics )

