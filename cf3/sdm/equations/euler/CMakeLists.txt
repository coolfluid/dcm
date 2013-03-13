# coolfluid_sdm_euler

list( APPEND coolfluid_sdm_equations_euler_files
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
  Terms1D.hpp
  Terms1D.cpp
  Terms2D.hpp
  Terms2D.cpp
  HLLE.hpp
  HLLE.cpp
  Roe.hpp
  Roe.cpp
  Rusanov.hpp
  Rusanov.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_equations_euler
                        SOURCES ${coolfluid_sdm_equations_euler_files}
                        LIBS    coolfluid_solver coolfluid_sdm_core coolfluid_physics_euler
                        PLUGINS Physics )

