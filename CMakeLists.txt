#Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> 
project(patchmatch)


# SETUP A MAKE CLEAN TO CLEAN EVERYTHING
#
#
set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
      "CMakeCache.txt;CMakeFiles;Makefile;cmake_install.cmake")

# DEFAULT BUILD: RELEASE
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()

if (NOT TARGET IIOLIB)
   include_directories("iio")
   add_subdirectory(iio)
endif() 
SET( LIBS ${LIBS} IIOLIB)

FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
   message("OPENMP FOUND")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

#
add_executable (patchmatch patchmatch_simple.cc img.cc)
target_link_libraries(patchmatch ${LIBS})

