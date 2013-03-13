list( APPEND coolfluid_sdm_tools_files
  LibTools.hpp
  LibTools.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_tools
                        SOURCES ${coolfluid_sdm_tools_files}
                        LIBS    coolfluid_solver_actions )


