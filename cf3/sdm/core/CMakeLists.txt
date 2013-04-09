list( APPEND coolfluid_sdm_core_files
  BC.hpp
  BC.cpp
  CellConnectivity.hpp
  CellConnectivity.cpp
  CombinedTermComputer.hpp
  CombinedTermComputer.cpp
  FaceConnectivity.hpp
  FaceConnectivity.cpp
  LegendreGaussLobatto.hpp
  LegendreGaussLobatto.cpp
  LibCore.hpp
  LibCore.cpp
  Metrics.hpp
  Metrics.cpp
  Reconstructions.hpp
  ShapeFunction.hpp
  ShapeFunction.cpp
  Tensorial.hpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_core
                        SOURCES ${coolfluid_sdm_core_files}
                        LIBS    coolfluid_mesh )

