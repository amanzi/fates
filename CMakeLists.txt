# Require cmake 3.1.3 or higher for c++11 support                                                                                                                                                           
cmake_minimum_required(VERSION 3.1.3)

set(GLOBAL LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(GLOBAL INCLUDE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)
set(GLOBAL RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
set(GLOBAL PROJECT_DIR ${CPIC_DIR})

# Source files for all binaries and libraries found under src                                                                                                                                               
add_subdirectory(shared)
#add_subdirectory(interface)
#add_subdirectory(biogeochem)
add_subdirectory(main)
#add_subdirectory(interface)
