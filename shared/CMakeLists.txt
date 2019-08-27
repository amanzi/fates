enable_language (C Fortran)
set(genf90_files shr_infnan_mod.F90.in shr_assert_mod.F90.in)

#process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR}
#  share_genf90_sources)

#sourcelist_to_parent(share_genf90_sources)

#list(APPEND share_sources "${share_genf90_sources}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)

set(ADDITIONAL_MAKE_CLEAN_FILES build)

get_filename_component (default_prefix "." ABSOLUTE)

include_directories(${NETCDF_DIR}
                    ${PIO_DIR})

set(share_sources
#  shr_infnan_mod.F90        
#  shr_assert_mod.F90
#  shr_file_mod.F90
#  shr_kind_mod.F90
#  shr_const_mod.F90
#  shr_sys_mod.F90
#  shr_log_mod.F90
#  shr_orb_mod.F90
#  shr_spfn_mod.F90
#  shr_strconvert_mod.F90
#  shr_cal_mod.F90
#  shr_nl_mod.F90
#  shr_precip_mod.F90
#  shr_string_mod.F90
#  shr_timer_mod.F90
#  shr_vmath_mod.F90
#  shr_wv_sat_mod.F90
  shr_kind_mod.F90
  shr_const_mod.F90
  shr_mpi_mod.F90
  shr_file_mod.F90
  shr_nl_mod.F90
  shr_infnan_mod.F90
  shr_string_mod.F90
  shr_sys_mod.F90
  shr_assert_mod.F90
  shr_log_mod.F90
  shr_timer_mod.F90
  shr_pio_mod.F90  
  abortutils.F90
  pio/fileutils.F90
  pio/ncdio_pio.F90
  pio/paramUtilMod.F90
)
configure_file(shr_assert.h shr_assert.h COPYONLY)
configure_file(mpif.h mpif.h COPYONLY)
configure_file(dtypes.h dtypes.h COPYONLY)

#include(add_betr_library)
#add_betr_library(betr_shr ${BETR_SHR_SOURCES})

# turn off all compiler warnings in clm code
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

#set(BETR_LIBRARIES betr_shr;${BETR_LIBRARIES} PARENT_SCOPE)
#set(BETR_LIBRARIES betr_shr;${BETR_LIBRARIES})


add_library(shared ${share_sources})
target_link_libraries(shared) 
install (TARGETS shared
       LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
       ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
