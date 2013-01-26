# coolfluid_sdm_advectiondiffusion

list( APPEND coolfluid_sdm_equations_advectiondiffusion_files
  LibAdvectionDiffusion.hpp
  LibAdvectionDiffusion.cpp
  AdvectionDiffusion1D.hpp
  AdvectionDiffusion1D.cpp
  Terms1D.hpp
  Terms1D.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_equations_advectiondiffusion
                        SOURCES ${coolfluid_sdm_equations_advectiondiffusion_files}
                        LIBS    coolfluid_solver coolfluid_sdm_core )

