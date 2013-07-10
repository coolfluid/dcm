add_subdirectory( core )
add_subdirectory( solver )
add_subdirectory( tools )
add_subdirectory( equations )

list( APPEND coolfluid_dcm_files
  LibDCM.hpp
  LibDCM.cpp
  Model.hpp
  Model.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm
                        SOURCES ${coolfluid_dcm_files}
                        LIBS    coolfluid_mesh_actions 
                                coolfluid_solver
                                coolfluid_solver_actions
                                coolfluid_physics
                                coolfluid_mesh_lagrangep2
                                coolfluid_dcm_core
                                coolfluid_dcm_tools )
