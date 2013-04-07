# sdm tests

#########################################################################

coolfluid_add_test( UTEST      utest-sdm-transformation
                    CPP        utest-sdm-transformation.cpp
                    LIBS       coolfluid_sdm )

coolfluid_add_test( UTEST      utest-sdm-faceconnectivity
                    CPP        utest-sdm-faceconnectivity.cpp
                    LIBS       coolfluid_sdm )

coolfluid_add_test( UTEST      utest-sdm-computerhs
                    CPP        utest-sdm-computerhs.cpp
                    LIBS       coolfluid_sdm coolfluid_mesh_gmsh
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-tools-computefieldgradientBR2
                    CPP        utest-sdm-tools-computefieldgradientBR2.cpp
                    LIBS       coolfluid_sdm coolfluid_mesh_gmsh
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-advectiondiffusion
                    CPP        utest-sdm-advectiondiffusion.cpp
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_advectiondiffusion coolfluid_sdm_solver_erkls
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-euler
                    CPP        utest-sdm-euler.cpp Term.hpp
                    PLUGINS    Physics
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_euler coolfluid_sdm_solver_erkls
                    MPI        1 )

coolfluid_add_test( UTEST      utest-sdm-lineuler-2d
                    CPP        utest-sdm-lineuler-2d.cpp
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_lineuler 
                    MPI        1
                    CONDITION  OFF )


coolfluid_add_test( UTEST      utest-sdm-linearadv1d
                    CPP        utest-sdm-linearadv1d.cpp
                    PLUGINS    Physics
                    LIBS       coolfluid_sdm_scalar coolfluid_mesh_gmsh coolfluid_mesh_tecplot coolfluid_physics_scalar
                    MPI        2 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-convection
                    CPP        utest-sdm-convection.cpp
                    LIBS       coolfluid_sdm coolfluid_mesh_gmsh coolfluid_mesh_tecplot coolfluid_physics_scalar
                    MPI        1 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-diffusion
                    CPP        utest-sdm-diffusion.cpp
                    LIBS       coolfluid_sdm coolfluid_mesh_gmsh coolfluid_sdm_equations_scalar
                    MPI        1 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-navierstokesmovingreference-2d
                    CPP        utest-sdm-navierstokesmovingreference-2d.cpp
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_navierstokesmovingreference 
                    CONDITION  OFF )

coolfluid_add_test( UTEST      utest-sdm-navierstokesmovingreference-3d
                    CPP        utest-sdm-navierstokesmovingreference-3d.cpp
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_navierstokesmovingreference 
                    CONDITION  OFF )

coolfluid_add_test( ATEST      atest-sdm-scalar-linadv-3d
                    PYTHON     atest-sdm-scalar-linadv-3d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_scalar
                    CONDITION  OFF )
                    
coolfluid_add_test( ATEST      atest-sdm-euler-shocktube-1d
                    PYTHON     atest-sdm-euler-shocktube-1d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_euler )

coolfluid_add_test( ATEST      atest-sdm-euler-shocktube-2d
                    PYTHON     atest-sdm-euler-shocktube-2d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_euler )

coolfluid_add_test( ATEST      atest-sdm-euler-cylinder-2d
                    PYTHON     atest-sdm-euler-cylinder-2d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_euler
                    MPI        2
                    DEPENDS    copy-resources )

coolfluid_add_test( ATEST      atest-sdm-lineuler-acousticpulse-2d
                    PYTHON     atest-sdm-lineuler-acousticpulse-2d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_lineuler
                    DEPENDS    copy-resources)

coolfluid_add_test( ATEST      atest-sdm-lineuler-wallbc-2d
                    PYTHON     atest-sdm-lineuler-wallbc-2d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_lineuler
                    DEPENDS    copy-resources)
                    
coolfluid_add_test( ATEST      atest-sdm-lineuler-acousticvorticityentropypulse-2d
                    PYTHON     atest-sdm-lineuler-acousticvorticityentropypulse-2d.py
                    LIBS       coolfluid_sdm_equations_lineuler
                    DEPENDS    copy-resources )
                    

coolfluid_add_test( ATEST      atest-sdm-lineuler-acousticpulse-3d
                    PYTHON     atest-sdm-lineuler-acousticpulse-3d.py
                    LIBS       coolfluid_sdm coolfluid_sdm_equations_lineuler
                    DEPENDS    copy-resources 
                    CONDITION  OFF )

coolfluid_add_test( ATEST      atest-sdm-lineuler-machcone-3d
                    PYTHON     atest-sdm-lineuler-machcone-3d.py
                    LIBS       coolfluid_sdm_equations_lineuler
                    DEPENDS    copy-resources
                    MPI        2 
                    CONDITION  OFF )

coolfluid_add_test( ATEST      atest-sdm-channelflow-movingref-2d
                    PYTHON     atest-sdm-channelflow-movingref-2d.py
                    LIBS       coolfluid_sdm_navierstokes coolfluid_sdm_equations_navierstokesmovingreference
                    MPI        1 
                    CONDITION  OFF )

#########################################################################
