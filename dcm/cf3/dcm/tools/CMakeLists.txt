list( APPEND coolfluid_dcm_tools_files
  LibTools.hpp
  LibTools.cpp
  ComputeFieldGradientBR2.hpp
  ComputeFieldGradientBR2.cpp
  CreateTermComputer.hpp
  CreateTermComputer.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_tools
                        SOURCES ${coolfluid_dcm_tools_files}
                        LIBS    coolfluid_solver_actions
                                coolfluid_dcm_core )


