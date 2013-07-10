list( APPEND coolfluid_dcm_core_files
  BC.hpp
  BC.cpp
  CellConnectivity.hpp
  CellConnectivity.cpp
  FaceConnectivity.hpp
  FaceConnectivity.cpp
  LegendreGaussEnd.hpp
  LegendreGaussEnd.cpp
  LibCore.hpp
  LibCore.cpp
  Metrics.hpp
  Metrics.cpp
  Reconstructions.hpp
  ShapeFunction.hpp
  ShapeFunction.cpp
  Tensorial.hpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_core
                        SOURCES ${coolfluid_dcm_core_files}
                        LIBS    coolfluid_mesh )

