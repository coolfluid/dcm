add_subdirectory( core )
add_subdirectory( solver )
add_subdirectory( tools )
add_subdirectory( equations )

list( APPEND coolfluid_sdm_files
  LibSDM.hpp
  LibSDM.cpp
  Model.hpp
  Model.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm
                        SOURCES ${coolfluid_sdm_files}
                        LIBS    coolfluid_mesh_actions 
                                coolfluid_solver
                                coolfluid_solver_actions
                                coolfluid_physics
                                coolfluid_mesh_lagrangep2
                                coolfluid_sdm_core )
