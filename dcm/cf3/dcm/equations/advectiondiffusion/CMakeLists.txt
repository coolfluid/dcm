# coolfluid_dcm_advectiondiffusion

list( APPEND coolfluid_dcm_equations_advectiondiffusion_files
  LibAdvectionDiffusion.hpp
  LibAdvectionDiffusion.cpp
  AdvectionDiffusion1D.hpp
  AdvectionDiffusion1D.cpp
  AdvectionDiffusion2D.hpp
  AdvectionDiffusion2D.cpp
  RightHandSide1D.hpp
  RightHandSide1D.cpp
  RightHandSide2D.hpp
  RightHandSide2D.cpp
  BCDirichlet.hpp
  BCDirichlet.cpp
  BCNeumann.hpp
  BCNeumann.cpp
  BCExtrapolation.hpp
  BCExtrapolation.cpp
  Terms.hpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_equations_advectiondiffusion
                        SOURCES ${coolfluid_dcm_equations_advectiondiffusion_files}
                        LIBS    coolfluid_solver coolfluid_dcm_core )

