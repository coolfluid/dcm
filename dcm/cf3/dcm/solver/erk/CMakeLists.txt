list( APPEND coolfluid_dcm_solver_erk_files 
  ButcherTableau.hpp
  ButcherTableau.cpp
  RungeKutta.hpp
  RungeKutta.cpp
  Types.hpp
  Types.cpp
  LibERK.hpp
  LibERK.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_solver_erk
                        SOURCES ${coolfluid_dcm_solver_erk_files}
                        LIBS    coolfluid_solver )
