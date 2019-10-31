enable_language (C Fortran)
set(genf90_files shr_infnan_mod.F90.in shr_assert_mod.F90.in)

#process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR}
#  share_genf90_sources)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)

set(ADDITIONAL_MAKE_CLEAN_FILES build)

get_filename_component (default_prefix "." ABSOLUTE)

include_directories(${NETCDF_INCLUDE_DIR}
                    ${PIO_INCLUDE_DIR})

set(parteh_sources
  PRTGenericMod.F90
  PRTAllometricCarbonMod.F90
  PRTLossFluxesMod.F90
)
#configure_file(mpif.h mpif.h COPYONLY)
#configure_file(dtypes.h dtypes.h COPYONLY)

#include(add_betr_library)
#add_betr_library(betr_shr ${BETR_SHR_SOURCES})

# turn off all compiler warnings in clm code
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

#set(BETR_LIBRARIES betr_shr;${BETR_LIBRARIES} PARENT_SCOPE)
#set(BETR_LIBRARIES betr_shr;${BETR_LIBRARIES})


add_library(parteh ${parteh_sources})
target_link_libraries(parteh) 
install (TARGETS parteh
       LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
       ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
