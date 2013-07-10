list( APPEND coolfluid_sdm_files
  LibSDM.hpp
  LibSDM.cpp
  BR2.hpp
  BR2.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm
                        SOURCES ${coolfluid_sdm_files}
                        LIBS    
                            coolfluid_dcm_equations_advectiondiffusion
                            coolfluid_dcm_equations_lineuler
                            coolfluid_dcm_equations_euler
                            coolfluid_dcm_equations_navierstokes
                            coolfluid_dcm_equations_les
                        PLUGINS DCM )

add_subdirectory( solver )