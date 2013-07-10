list( APPEND coolfluid_sdm_solver_optim_erkls_files 
  LibOptimERKLS.hpp
  LibOptimERKLS.cpp
  Types.hpp
  Types.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_solver_optim_erkls
                        SOURCES ${coolfluid_sdm_solver_optim_erkls_files}
                        LIBS    coolfluid_dcm_solver_erkls
                        PLUGINS DCM )
