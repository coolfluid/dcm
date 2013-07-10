list( APPEND coolfluid-machcone_files  coolfluid-machcone.cpp  )
list( APPEND coolfluid-machcone_cflibs ${CF3_KERNEL_LIBS} coolfluid_sdm_lineuler )
set( coolfluid-machcone_condition OFF )
coolfluid_add_application( coolfluid-machcone )
