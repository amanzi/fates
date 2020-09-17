project (FATES)
enable_language (C Fortran)
# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)

set(ADDITIONAL_MAKE_CLEAN_FILES build)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)

MESSAGE( STATUS "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX:     " ${CMAKE_INSTALL_PREFIX} )
MESSAGE( STATUS "LIBRARY_OUTPUT_DIRECTORY:     " ${LIBRARY_OUTPUT_DIRECTORY} )
MESSAGE( STATUS "RUNTIME_OUTPUT_DIRECTORY:     " ${RUNTIME_OUTPUT_DIRECTORY} )
MESSAGE( STATUS "INCLUDE_OUTPUT_DIRECTORY:     " ${INCLUDE_OUTPUT_DIRECTORY} )
MESSAGE( STATUS "CMAKE_Fortran_MODULE_DIRECTORY:     " ${CMAKE_Fortran_MODULE_DIRECTORY} )


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

include_directories(${SHARED_SOURCE_DIR}
                    ${CMAKE_INSTALL_PREFIX}/include )

                  set(fates_src_files
                    FatesConstantsMod.F90
                    FatesGlobals.F90
                    FatesUtilsMod.F90
                    FatesIODimensionsMod.F90
                    FatesIOVariableKindMod.F90
                    FatesHistoryInterfaceMod.F90
                    FatesHistoryVariableType.F90
                    FatesHydraulicsMemMod.F90
                    FatesIntegratorsMod.F90
                    FatesInterfaceMod.F90
                    FatesInventoryInitMod.F90
                    FatesParameterDerivedMod.F90
                    FatesParametersInterface.F90
                    FatesRestartInterfaceMod.F90
                    FatesRestartVariableType.F90
                    FatesSizeAgeTypeIndicesMod.F90
                    FatesSynchronizedParamsMod.F90
                    ChecksBalancesMod.F90
                    EDInitMod.F90
                    EDMainMod.F90
                    EDParamsMod.F90
                    EDPftvarcon.F90
                    EDTypesMod.F90
                    histFileMod.F90
                    ../parteh/PRTGenericMod.F90
                    ../parteh/PRTAllometricCarbonMod.F90
                    ../parteh/PRTLossFluxesMod.F90
                    ../fire/SFMainMod.F90
                    ../fire/SFParamsMod.F90
                    ../biogeochem/EDCanopyStructureMod.F90
                    ../biogeochem/EDCohortDynamicsMod.F90
                    ../biogeochem/EDLoggingMortalityMod.F90
                    ../biogeochem/EDMortalityFunctionsMod.F90
                    ../biogeochem/EDPatchDynamicsMod.F90
                    ../biogeochem/EDPhysiologyMod.F90
                    ../biogeochem/FatesAllometryMod.F90
                    ../biogeochem/FatesLitterMod.F90
                    ../biogeophys/EDAccumulateFluxesMod.F90
                    ../biogeophys/EDBtranMod.F90
                    ../biogeophys/FatesPlantHydraulicsMod.F90
                    ../biogeophys/EDSurfaceAlbedoMod.F90
                    ../biogeophys/FatesPlantRespPhotosynthMod.F90
                    ../biogeophys/FatesBstressMod.F90
                    ../interface/atsfates_interfaceMod.F90
                    ../interface/atsfates_paraminterfaceMod.F90
                    )

add_library(fates ${fates_src_files})
target_link_libraries(fates) 
install (TARGETS fates
       LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
       ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

