# coolfluid_sdm_scalar

list( APPEND coolfluid_sdm_equations_scalar_files 
  LibScalar.hpp
  LibScalar.cpp 
)
coolfluid_mark_not_orphan(
# scalar 1D
  LinearAdvection1D.hpp
  LinearAdvection1D.cpp
  Diffusion1D.hpp
  Diffusion1D.cpp
# scalar 2D
  LinearAdvection2D.hpp
  LinearAdvection2D.cpp
  Diffusion2D.hpp
  Diffusion2D.cpp
  RotationAdvection2D.hpp
  RotationAdvection2D.cpp
# scalar 3D
  LinearAdvection3D.hpp
  LinearAdvection3D.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_equations_scalar
                        SOURCES ${coolfluid_sdm_equations_scalar_files}
                        LIBS    coolfluid_sdm_core )

