coolfluid3_add_library( TARGET    coolfluid_sdm_br2
                        SOURCES
                            LibBR2.hpp
                            LibBR2.cpp
                            BR2.hpp
                            BR2.cpp
                        LIBS    
                            coolfluid_dcm_equations_advectiondiffusion
                            coolfluid_dcm_equations_lineuler
                            coolfluid_dcm_equations_euler
                            coolfluid_dcm_equations_navierstokes
                            coolfluid_dcm_equations_les
                        PLUGINS DCM )
