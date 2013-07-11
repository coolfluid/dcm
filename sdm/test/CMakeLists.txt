# sdm tests

#########################################################################

coolfluid_add_test( UTEST      utest-dcm-transformation
                    CPP        utest-dcm-transformation.cpp
                    LIBS       coolfluid_dcm )

coolfluid_add_test( UTEST      utest-dcm-faceconnectivity
                    CPP        utest-dcm-faceconnectivity.cpp
                    LIBS       coolfluid_dcm )

coolfluid_add_test( UTEST      utest-dcm-tools-computefieldgradientBR2
                    CPP        utest-dcm-tools-computefieldgradientBR2.cpp
                    LIBS       coolfluid_dcm coolfluid_mesh_gmsh coolfluid_mesh
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-computerhs
                    CPP        utest-sdm-computerhs.cpp Term.hpp
                    LIBS       coolfluid_dcm coolfluid_sdm_br2 coolfluid_mesh_gmsh
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-advectiondiffusion
                    CPP        utest-sdm-advectiondiffusion.cpp
                    LIBS       coolfluid_dcm coolfluid_sdm_br2 coolfluid_dcm_equations_advectiondiffusion coolfluid_dcm_solver_erkls
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-euler
                    CPP        utest-sdm-euler.cpp Term.hpp
                    PLUGINS    Physics
                    LIBS       coolfluid_dcm coolfluid_sdm_br2 coolfluid_dcm_equations_euler coolfluid_dcm_solver_erkls
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-navierstokes
                    CPP        utest-sdm-navierstokes.cpp Term.hpp
                    PLUGINS    Physics
                    LIBS       coolfluid_dcm coolfluid_sdm_br2 coolfluid_dcm_equations_navierstokes coolfluid_dcm_solver_erkls
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-les
                    CPP        utest-sdm-les.cpp
                    PLUGINS    Physics
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_les
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-lineuler-2d
                    CPP        utest-sdm-lineuler-2d.cpp
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_lineuler
                    MPI        1
                    CONDITION  OFF )


coolfluid_add_test( UTEST      utest-sdm-linearadv1d
                    CPP        utest-sdm-linearadv1d.cpp
                    PLUGINS    Physics
                    LIBS       coolfluid_mesh_gmsh coolfluid_mesh_tecplot coolfluid_physics_scalar
                    MPI        2 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-convection
                    CPP        utest-sdm-convection.cpp
                    LIBS       coolfluid_sdm_br2 coolfluid_mesh_gmsh coolfluid_mesh_tecplot coolfluid_physics_scalar
                    MPI        1 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-diffusion
                    CPP        utest-sdm-diffusion.cpp
                    LIBS       coolfluid_sdm_br2 coolfluid_mesh_gmsh coolfluid_dcm_equations_scalar
                    MPI        1 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-navierstokesmovingreference-2d
                    CPP        utest-sdm-navierstokesmovingreference-2d.cpp
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_navierstokesmovingreference
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-navierstokesmovingreference-3d
                    CPP        utest-sdm-navierstokesmovingreference-3d.cpp
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_navierstokesmovingreference
                    CONDITION  OFF )

coolfluid_add_test( ATEST      atest-sdm-scalar-linadv-3d
                    PYTHON     atest-sdm-scalar-linadv-3d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_scalar
                    CONDITION  OFF )
                    
coolfluid_add_test( ATEST      atest-sdm-euler-shocktube-1d
                    PYTHON     atest-sdm-euler-shocktube-1d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_euler )

coolfluid_add_test( ATEST      atest-sdm-euler-shocktube-2d
                    PYTHON     atest-sdm-euler-shocktube-2d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_euler )

coolfluid_add_test( ATEST      atest-sdm-euler-cylinder-2d
                    PYTHON     atest-sdm-euler-cylinder-2d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_euler
                    MPI        2
                    DEPENDS    copy-resources )

coolfluid_add_test( ATEST      atest-sdm-lineuler-acousticpulse-2d
                    PYTHON     atest-sdm-lineuler-acousticpulse-2d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_lineuler
                    DEPENDS    copy-resources)

coolfluid_add_test( ATEST      atest-sdm-lineuler-wallbc-2d
                    PYTHON     atest-sdm-lineuler-wallbc-2d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_lineuler
                    DEPENDS    copy-resources)
                    
coolfluid_add_test( ATEST      atest-sdm-lineuler-acousticvorticityentropypulse-2d
                    PYTHON     atest-sdm-lineuler-acousticvorticityentropypulse-2d.py
                    LIBS       coolfluid_dcm_equations_lineuler coolfluid_sdm_solver_optim_erkls
                    DEPENDS    copy-resources )
                    

coolfluid_add_test( ATEST      atest-sdm-lineuler-acousticpulse-3d
                    PYTHON     atest-sdm-lineuler-acousticpulse-3d.py
                    LIBS       coolfluid_sdm_br2 coolfluid_dcm_equations_lineuler
                    DEPENDS    copy-resources 
                    CONDITION  OFF )

coolfluid_add_test( ATEST      atest-sdm-lineuler-machcone-3d
                    PYTHON     atest-sdm-lineuler-machcone-3d.py
                    LIBS       coolfluid_dcm_equations_lineuler
                    DEPENDS    copy-resources
                    MPI        2 
                    CONDITION  OFF )

coolfluid_add_test( ATEST      atest-sdm-channelflow-movingref-2d
                    PYTHON     atest-sdm-channelflow-movingref-2d.py
                    LIBS       coolfluid_dcm_equations_navierstokesmovingreference
                    MPI        1 
                    CONDITION  OFF )

#########################################################################
