##############################################################################
# Discontinuous Collocation Methods plugin
##############################################################################

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )

project( discontinuous_collocation_methods CXX )

coolfluid_define_plugin( DCM ${CMAKE_CURRENT_SOURCE_DIR} )

if( CF3_PLUGIN_DCM )

  coolfluid_find_orphan_files()

  include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

  add_subdirectory( cf3 )   # library sources
    
  # add_subdirectory( test )  # testing tree

  # add_subdirectory( executables ) # executables

endif()
