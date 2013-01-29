# coolfluid_sdm_navierstokes

list( APPEND coolfluid_sdm_equations_navierstokes_files
  BCWall.hpp
  BCWall.cpp
  LibNavierStokes.hpp
  LibNavierStokes.cpp
  NavierStokes1D.hpp
  NavierStokes1D.cpp
  NavierStokes2D.hpp
  NavierStokes2D.cpp
  Terms1D.hpp
  Terms1D.cpp
  Terms2D.hpp
  Terms2D.cpp
  HLLE.cpp
  HLLE.hpp
  Roe.cpp
  Roe.hpp
  Rusanov.cpp
  Rusanov.hpp
)
coolfluid_mark_not_orphan(
  BCSubsonicInlet2D.hpp
  BCSubsonicInlet2D.cpp
  BCSubsonicOutlet2D.hpp
  BCSubsonicOutlet2D.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_equations_navierstokes
                        SOURCES ${coolfluid_sdm_equations_navierstokes_files}
                        LIBS    coolfluid_sdm_equations_euler coolfluid_physics_navierstokes
                        PLUGINS Physics )

