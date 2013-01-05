list( APPEND coolfluid_sdm_solver_erk_files 
  ButcherTableau.hpp
  ButcherTableau.cpp
  RungeKutta.hpp
  RungeKutta.cpp
  Types.hpp
  Types.cpp
  LibERK.hpp
  LibERK.cpp
)

coolfluid3_add_library( TARGET  coolfluid_sdm_solver_erk
                        SOURCES ${coolfluid_sdm_solver_erk_files}
                        LIBS    coolfluid_solver )
