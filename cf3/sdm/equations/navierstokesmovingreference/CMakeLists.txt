# coolfluid_sdm_navierstokesmovingreference

list( APPEND coolfluid_sdm_equations_navierstokesmovingreference_files
  LibNavierStokesMovingReference.hpp
  LibNavierStokesMovingReference.cpp 
)
coolfluid_mark_not_orphan (
  BCSubsonicInlet2D.hpp
  BCSubsonicInlet2D.cpp
  BCSubsonicInlet3D.hpp
  BCSubsonicInlet3D.cpp
  BCSubsonicOutlet2D.hpp
  BCSubsonicOutlet2D.cpp
  BCSubsonicOutlet3D.hpp
  BCSubsonicOutlet3D.cpp
  Convection2D.hpp
  Convection2D.cpp
  Convection3D.hpp
  Convection3D.cpp
  Source2D.hpp
  Source2D.cpp
  Source3D.hpp
  Source3D.cpp
  Diffusion2D.hpp
  Diffusion2D.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_equations_navierstokesmovingreference
                        SOURCES ${coolfluid_sdm_equations_navierstokesmovingreference_files}
                        LIBS    coolfluid_sdm_core coolfluid_physics_navierstokes
                        PLUGINS Physics )

