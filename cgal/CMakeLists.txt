project( preprocessing )

cmake_minimum_required(VERSION 2.8.11)

#IF(WIN32)
#IF(MSVC)
#ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
#ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
#ENDIF(MSVC)
#ELSE(WIN32)
#SET(CMAKE_CXX_FLAGS "-W -Wall -Wno-unused")
#ENDIF(WIN32)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  

endif()

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)

if (NOT EIGEN3_FOUND)
  find_package(LAPACK)
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
  endif(LAPACK_FOUND)
else()
  include( ${EIGEN3_USE_FILE} )
endif()

# include for local directory

# include for local package


# Creating entries for target: executable
# ############################

add_executable( preprocessing preprocessing.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS preprocessing )

# Link the executable to CGAL and third-party libraries
target_link_libraries(preprocessing ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

