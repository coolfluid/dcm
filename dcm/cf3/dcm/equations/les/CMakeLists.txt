# coolfluid_dcm_les

list( APPEND coolfluid_dcm_equations_les_files
  LibLES.hpp
  LibLES.cpp
  LES2D.hpp
  LES2D.cpp
  EddyViscosityModel.hpp
  EddyViscosityModel.cpp
  RightHandSide2D.hpp
  RightHandSide2D.cpp
  Terms.hpp
  ComputeSubFilterScale.hpp
  ComputeSubFilterScale.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_equations_les
                        SOURCES ${coolfluid_dcm_equations_les_files}
                        LIBS    coolfluid_dcm_equations_navierstokes
                        PLUGINS Physics )

