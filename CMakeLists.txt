cmake_minimum_required(VERSION 3.16)
project(ULocator VERSION 0.3.0 LANGUAGES C CXX)
enable_testing()

# Some options
option(WRAP_PYTHON "Compile the Python bindings" OFF)
option(BUILD_UUSS "Compile the UUSS examples" OFF)
option(USE_UMPS "Compile with the UMPS logging library" ON)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(GTest REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Eigen3 3.3 REQUIRED Dense)
find_package(GeographicLib REQUIRED)
#find_package(TTimes REQUIRED)
find_package(HDF5 COMPONENTS C REQUIRED)
#find_package(MCMC REQUIRED)
find_package(NLopt REQUIRED)
find_package(TBB)
set(ONEDAL_USE_DPCPP no) 
set(ONEDAL_INTERFACE yes)
find_package(oneDAL REQUIRED)
if (BUILD_UUSS)
   add_compile_definitions(WITH_HDF5)
   find_package(HDF5 REQUIRED)
   find_package(nlohmann_json 3.2.0 REQUIRED)
endif()
if (${TBB_FOUND})
   message("Using TBB")
   add_compile_definitions(USE_TBB)
endif()
if (${USE_UMPS})
   message("Will compile with UMPS logging")
   add_compile_definitions(WITH_UMPS)
endif()

# Versioning information
configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/uLocator/version.hpp)

set(LIBRARY_SRC src/version.cpp
                src/homogeneous.cpp
                src/arrival.cpp
                #src/direct.cpp
                #src/directOptions.cpp
                src/firstArrivalRayTracer.cpp
                #src/hamiltonianMonteCarlo.cpp
                #src/hamiltonianMonteCarloOptions.cpp
                #src/lineSearchOptions.cpp
                src/nlopt.cpp
                src/nloptOptions.cpp
                src/origin.cpp
                src/quarry.cpp
                src/staticCorrection.cpp
                src/sourceSpecificStationCorrection.cpp
                src/station.cpp
                src/topography.cpp
                src/travelTimeCalculatorMap.cpp
                src/uussTravelTimeCaclulator.cpp
                src/objectiveFunctions/l1.cpp
                src/objectiveFunctions/lp.cpp
                src/objectiveFunctions/doubleDifferenceL1.cpp
                src/objectiveFunctions/leastSquares.cpp
                src/objectiveFunctions/objectiveFunction.cpp
                src/position/wgs84.cpp)
if (NOT ${USE_UMPS})
   set(LIBRARY_SRC ${LIBRARY_SRC} src/logging/standardOut.cpp)
endif()
if (${BUILD_UUSS})
#    set(UUSS_SRC src/topography.cpp)
endif()
set(BUILD_SHARED_LIBS YES)
add_library(uLocator ${LIBRARY_SRC} ${UUSS_SRC})

set_target_properties(uLocator PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO) 
target_link_libraries(uLocator
                      PRIVATE Threads::Threads spdlog::spdlog Eigen3::Eigen
                              hdf5::hdf5 NLopt::nlopt TBB::tbb hdf5::hdf5
                              Boost::program_options 
                              #${MCMC_LIBRARY} ttimes
                              GeographicLib::GeographicLib ${UMPS_LIBRARY}
                              ${EIKONALXX_LIBRARY} ${DAL_LIBS})
target_include_directories(uLocator
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE ${UMPS_INCLUDE_DIR} ${MCMC_INCLUDE_DIR}
                                   ${NLOPT_INCLUDE_DIRS} ${EIKONALXX_INCLUDE_DIR}
                                   ${oneDAL_INCLUDE_DIRS} hdf5::hdf5 NLopt::nlopt
                                   #ttimes
                           PUBLIC $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include/umps>)


##########################################################################################
#                                         Unit Tests                                     #
##########################################################################################
#file(COPY ${CMAKE_SOURCE_DIR}/testing/data DESTINATION .)
set(TEST_SRC
    testing/main.cpp
    testing/arrival.cpp
    testing/position.cpp
    testing/topography.cpp
    testing/nlopt.cpp
    testing/homogeneous.cpp
    testing/example.cpp
    )
add_executable(unitTests ${TEST_SRC})
set_target_properties(unitTests PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_link_libraries(unitTests PRIVATE uLocator ${GTEST_BOTH_LIBRARIES} ${UMPS_LIBRARY} Threads::Threads)
target_include_directories(unitTests
                           PRIVATE ${GTEST_INCLUDE_DIRS}
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                           PRIVATE $<BUILD_INTERFACE:${UMPS_INCLUDE_DIR}>)
add_test(NAME unitTests 
         COMMAND unitTests)

##########################################################################################
#                                        Executables                                     #
##########################################################################################
set(UTAH_SRC)
add_executable(ulocate src/ulocate.cpp)
set_target_properties(ulocate PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO) 
target_link_libraries(ulocate
                      PRIVATE uLocator Threads::Threads
                      nlohmann_json::nlohmann_json Boost::program_options TBB::tbb
                      ${UMPS_LIBRARY})
target_include_directories(ulocate
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE $<BUILD_INTERFACE:${UMPS_INCLUDE_DIR}>)

add_executable(computeCorrections src/computeCorrections.cpp)
set_target_properties(computeCorrections PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_include_directories(computeCorrections
                           PRIVATE Boost::program_options ${oneDAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
message(" ${oneDAL_INCLUDE_DIRS}")
message(" ${DAL_LIBS}")
target_link_libraries(computeCorrections
                      ${DAL_LIBS} Boost::program_options uLocator)

add_executable(estimateVelocities
               src/estimateVelocities.cpp)
set_target_properties(estimateVelocities PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_link_libraries(estimateVelocities
                      uLocator Boost::program_options TBB::tbb ${UMPS_LIBRARY}
                      NLopt::nlopt)
target_include_directories(estimateVelocities
                           PRIVATE Boost::program_options ${CMAKE_SOURCE_DIR}/include
                                   ${UMPS_INCLUDE_DIR})

add_executable(evaluateCorrections
               src/evaluateCorrections.cpp)
set_target_properties(evaluateCorrections PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_link_libraries(evaluateCorrections uLocator)
target_include_directories(evaluateCorrections
                           PRIVATE Boost::program_options ${CMAKE_SOURCE_DIR}/include
                                   ${UMPS_INCLUDE_DIR})


if (${WRAP_PYTHON})
   message("Will build Python bindings")
   find_package(pybind11 REQUIRED)
endif()

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "ben.baker@utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A Hamiltonian Monte Carlo earthquake location utility.")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_COMPONENTS_ALL libraries headers)
if (WIN32)
    set(CPACK_GENERATOR ZIP WIX)
elseif (APPLE)
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CPACK_GENERATOR TGZ RPM)
else()
    set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_SOURCE_IGNORE_FILES
  /\\.git/
  \\.swp
  \\.orig
  /CMakeLists\\.txt\\.user
  /private/
)
include(CPack) # Put this last!
