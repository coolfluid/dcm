list( APPEND coolfluid_sdm_core_files
  BC.hpp
  BC.cpp
  CellConnectivity.hpp
  CellConnectivity.cpp
  CombinedTermComputer.hpp
  CombinedTermComputer.cpp
  FaceConnectivity.hpp
  FaceConnectivity.cpp
  LagrangeLocally1D.hpp
  LagrangeLocally1D.cpp
  LibCore.hpp
  LibCore.cpp
  Metrics.hpp
  Metrics.cpp
  Reconstructions.hpp
  ShapeFunction.hpp
  ShapeFunction.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_core
                        SOURCES ${coolfluid_sdm_core_files}
                        LIBS    coolfluid_mesh )

