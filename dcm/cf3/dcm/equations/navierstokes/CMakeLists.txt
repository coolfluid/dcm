# coolfluid_dcm_navierstokes

list( APPEND coolfluid_dcm_equations_navierstokes_files
  BCWall.hpp
  BCWall.cpp
  ComputePrimitiveVariables.hpp
  ComputePrimitiveVariables.cpp
  LibNavierStokes.hpp
  LibNavierStokes.cpp
  NavierStokes1D.hpp
  NavierStokes1D.cpp
  NavierStokes2D.hpp
  NavierStokes2D.cpp
  RightHandSide1D.hpp
  RightHandSide1D.cpp
  RightHandSide2D.hpp
  RightHandSide2D.cpp
  InitBlasius.hpp
  InitBlasius.cpp
  Terms.hpp
)
coolfluid_mark_not_orphan(
  BCSubsonicInlet2D.hpp
  BCSubsonicInlet2D.cpp
  BCSubsonicOutlet2D.hpp
  BCSubsonicOutlet2D.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_equations_navierstokes
                        SOURCES ${coolfluid_dcm_equations_navierstokes_files}
                        LIBS    coolfluid_dcm_equations_euler coolfluid_physics_navierstokes
                        PLUGINS Physics )

