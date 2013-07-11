list( APPEND coolfluid_dcm_solver_erkls_files 
  ThreeSstar.hpp
  ThreeSstar.cpp
  ThreeSstarCoeffs.hpp
  ThreeSstarCoeffs.cpp
  TwoSstar.hpp
  TwoSstar.cpp
  TwoSstarCoeffs.hpp
  TwoSstarCoeffs.cpp
  LibERKLS.hpp
  LibERKLS.cpp
)

coolfluid3_add_library( TARGET  coolfluid_dcm_solver_erkls
                        SOURCES ${coolfluid_dcm_solver_erkls_files}
                        LIBS    coolfluid_solver )
