cmake_minimum_required(VERSION 3.16)
project(uLocator VERSION 0.4.0 LANGUAGES C CXX)
enable_testing()

# Some options
option(WRAP_PYTHON "Compile the Python bindings" OFF)
option(BUILD_UUSS "Compile the UUSS examples" OFF)
option(USE_UMPS "Compile with the UMPS logging library" ON)
option(USE_EIKONALXX "Link against the ray tracer in the eikonal-solver software" ON)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(GTest REQUIRED)
#find_package(Catch2)
find_package(GeographicLib REQUIRED)

find_package(prima)
if (${BUILD_UUSS})
    find_package(pagmo REQUIRED)
else()
    find_package(pagmo)
endif()

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Threads REQUIRED)
find_package(NLopt REQUIRED)
find_package(TBB)
#find_package(zfp)
set(ONEDAL_LINK dynamic)
set(ONEDAL_USE_DPCPP no) 
set(ONEDAL_INTERFACE yes)
set(ONEDAL_SET_TBB_MANUALLY yes)
find_package(oneDAL REQUIRED)
if (BUILD_UUSS)
   find_package(HDF5 COMPONENTS C REQUIRED)
   find_package(nlohmann_json 3.2.0 REQUIRED)
else()
   find_package(HDF5 COMPONENTS C)
endif()
if (${HDF5_FOUND})
   message("HDF5 found")
   add_compile_definitions(WITH_HDF5)
endif()
if (${TBB_FOUND})
   message("Using TBB")
   add_compile_definitions(USE_TBB)
endif()
if (${USE_UMPS})
   message("Will compile with UMPS logging")
   add_compile_definitions(WITH_UMPS)
else()
   find_package(spdlog REQUIRED)
endif()
if (${prima_FOUND})
   message("Found prima")
   add_compile_definitions(WITH_PRIMA)
endif()
if (${pagmo_FOUND})
   message("Found pagmo")
   add_compile_definitions(WITH_PAGMO)
endif()
if (${NLopt_FOUND})
   message("NLOpt found")
   add_compile_definitions(WITH_NLOPT)
endif()
if (${USE_EIKONALXX})
   message("Will compile with EikonalXX ray tracer")
   add_compile_definitions(WITH_EIKONALXX)
else()
   message("Will compile internal ray tracer")
endif()

# Versioning information
configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/uLocator/version.hpp)

set(LIBRARY_SRC lib/version.cpp
                lib/arrival.cpp
                lib/origin.cpp
                lib/optimizers/optimizer.cpp
                lib/optimizers/originTime.cpp
                lib/rayTracer.cpp
                lib/station.cpp
                lib/travelTimeCalculator.cpp
                lib/travelTimeCalculatorMap.cpp
                lib/topography/constant.cpp
                lib/topography/gridded.cpp
                lib/position/geographicPoint.cpp
                lib/position/wgs84.cpp
                lib/corrections/static.cpp
                lib/corrections/sourceSpecific.cpp
                )
set(UTAH_SRC lib/uussRayTracer.cpp
             lib/position/knownUtahQuarry.cpp
             lib/position/knownUtahEvent.cpp
             lib/position/knownYNPEvent.cpp
             lib/position/utahRegion.cpp
             lib/position/utahQuarry.cpp
             lib/position/ynpRegion.cpp)
set(LIBRARY_SRC ${LIBRARY_SRC} ${UTAH_SRC})
if (${pagmo_FOUND})
   set(LIBRARY_SRC ${LIBRARY_SRC}
       lib/optimizers/pagmo/particleSwarm.cpp)
endif()
if (${prima_FOUND})
   set(LIBRARY_SRC ${LIBRARY_SRC}
      # src/optimizers/prima/cobyla.cpp
      lib/optimizers/prima/boundOptimizationByQuadraticApproximation.cpp
      )
endif()
if (${NLopt_FOUND})
   set(LIBRARY_SRC ${LIBRARY_SRC}
       lib/optimizers/nlopt/boundOptimizationByQuadraticApproximation.cpp
       lib/optimizers/nlopt/dividedRectangles.cpp
       lib/optimizers/nlopt/methodOfMovingAsymptotes.cpp
       lib/optimizers/nlopt/stochasticGradientOptimization.cpp)
endif()
if (NOT ${USE_UMPS})
   set(LIBRARY_SRC ${LIBRARY_SRC} lib/logging/standardOut.cpp)
endif()
if (NOT ${USE_EIKONALXX})
   set(LIBRARY_SRC ${LIBRARY_SRC}
       lib/ray/layerSolver.cpp
       lib/ray/point2d.cpp
       lib/ray/segment2d.cpp
       lib/ray/path2d.cpp)
endif()
set(BUILD_SHARED_LIBS YES)
add_library(uLocator ${LIBRARY_SRC} ${UUSS_SRC})

set_target_properties(uLocator PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO) 
target_link_libraries(uLocator
                      PRIVATE ${HDF5_C_LIBRARIES}
                              NLopt::nlopt
                              TBB::tbb
                              Boost::program_options Boost::boost
                              GeographicLib::GeographicLib
                              oneDAL::onedal oneDAL::onedal_core oneDAL::onedal_parameters oneDAL::onedal_thread
                              ${UMPS_LIBRARY} ${EIKONALXX_LIBRARY})
target_include_directories(uLocator
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE ${CMAKE_SOURCE_DIR}/lib
                                   ${UMPS_INCLUDE_DIR} ${MCMC_INCLUDE_DIR}
                                   ${NLOPT_INCLUDE_DIRS} ${EIKONALXX_INCLUDE_DIR}
                                   ${oneDAL_INCLUDE_DIRS} ${HDF5_C_INCLUDE_DIRS}
                           PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/umps>
                           PUBLIC $<INSTALL_INTERFACE:include>)

if (${pagmo_FOUND})
   target_include_directories(uLocator PRIVATE Pagmo::pagmo)
   target_link_libraries(uLocator PRIVATE Pagmo::pagmo)
endif()
if (${prima_FOUND})
   target_include_directories(uLocator PRIVATE prima::primac)
   target_link_libraries(uLocator PRIVATE prima::primac)
endif()
if (${NLopt_FOUND})
   target_include_directories(uLocator PRIVATE NLopt::nlopt)
   target_link_libraries(uLocator PRIVATE NLopt::nlopt)
endif() 

##########################################################################################
#                                         Unit Tests                                     #
##########################################################################################
#file(COPY ${CMAKE_SOURCE_DIR}/testing/data DESTINATION .)
set(TEST_SRC
    testing/main.cpp
    testing/position.cpp
    testing/station.cpp
    testing/arrival.cpp
    testing/origin.cpp
    testing/topography.cpp
    testing/rayTracer.cpp
    testing/corrections.cpp
    testing/objectiveFunctions.cpp
    testing/originTime.cpp
    testing/travelTimeDatabase.cpp
    #testing/halfSpace.cpp
    #testing/nlopt.cpp
    #testing/homogeneous.cpp
    #testing/example.cpp
    )
if (${pagmo_FOUND})
   set(TEST_SRC ${TEST_SRC} testing/pagmo.cpp)
endif()
if (${NLopt_FOUND})
   set(TEST_SRC ${TEST_SRC} testing/nlopt.cpp)
endif()

add_executable(unitTests ${TEST_SRC})
set_target_properties(unitTests PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_link_libraries(unitTests PRIVATE uLocator ${GTEST_BOTH_LIBRARIES} ${UMPS_LIBRARY}
                      oneDAL::onedal oneDAL::onedal_core oneDAL::onedal_parameters oneDAL::onedal_thread Boost::boost Threads::Threads)# nlohmann_json::nlohmann_json)
target_include_directories(unitTests
                           PRIVATE ${GTEST_INCLUDE_DIRS}
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib>
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                           PRIVATE $<BUILD_INTERFACE:${UMPS_INCLUDE_DIR}>)
add_test(NAME unitTests 
         COMMAND unitTests)

##########################################################################################
#                                        Executables                                     #
##########################################################################################
if (${BUILD_UUSS})
   add_executable(ulocate src/ulocate.cpp)
   set_target_properties(ulocate PROPERTIES
                         CXX_STANDARD 20
                         CXX_STANDARD_REQUIRED YES 
                         CXX_EXTENSIONS NO) 
   target_link_libraries(ulocate
                         PRIVATE uLocator Threads::Threads
                         Boost::program_options Boost::boost TBB::tbb
                         ${UMPS_LIBRARY})
   target_include_directories(ulocate
                              PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                              PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib>
                              PRIVATE $<BUILD_INTERFACE:${UMPS_INCLUDE_DIR}>)
endif()

# DOWNLOAD UTAH/YNP topography:
#file(DOWNLOAD "https://docs.google.com/uc?export=download&id=1WdZ5CWsCkWr19lXlQd7S_fLuuW_SCUGU" utahTopography.h5)
#file(DOWNLOAD "https://docs.google.com/uc?export=download&id=1HjwoDcWPz4Pj1JwJTijRvXH85hevMs_O" ynpTopography.h5)
if (${BUILD_UUSS})
   add_executable(computeCorrections src/computeCorrections.cpp)
   set_target_properties(computeCorrections PROPERTIES
                         CXX_STANDARD 20
                         CXX_STANDARD_REQUIRED YES
                         CXX_EXTENSIONS NO)
   target_include_directories(computeCorrections
                              PRIVATE Boost::program_options ${oneDAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/lib)
   target_link_libraries(computeCorrections
                         ${DAL_LIBS} Boost::program_options uLocator)
endif()

if (${BUILD_UUSS})
   add_executable(estimateVelocities src/estimateVelocities.cpp)
   set_target_properties(estimateVelocities PROPERTIES
                         CXX_STANDARD 20
                         CXX_STANDARD_REQUIRED YES
                         CXX_EXTENSIONS NO)
   target_link_libraries(estimateVelocities
                         uLocator Boost::program_options TBB::tbb ${UMPS_LIBRARY}
                         NLopt::nlopt prima::primac)
   target_include_directories(estimateVelocities
                           PRIVATE Boost::program_options
                                   ${CMAKE_SOURCE_DIR}/include
                                   ${CMAKE_SOURCE_DIR}/lib
                                   ${UMPS_INCLUDE_DIR} NLopt::nlopt prima::primac)
endif()
#add_executable(evaluateCorrections
#               src/evaluateCorrections.cpp)
#set_target_properties(evaluateCorrections PROPERTIES
#                      CXX_STANDARD 20
#                      CXX_STANDARD_REQUIRED YES
#                      CXX_EXTENSIONS NO)
#target_link_libraries(evaluateCorrections uLocator)
#target_include_directories(evaluateCorrections
#                           PRIVATE Boost::program_options ${CMAKE_SOURCE_DIR}/include
#                                   ${UMPS_INCLUDE_DIR})

if (${WRAP_PYTHON})
   message("Will build Python bindings")
   find_package(pybind11 REQUIRED)
   add_library(pyulocator MODULE
               python/pyulocator.cpp
               python/dataStructures.cpp
               python/position.cpp
               python/topography.cpp
               python/rayTracers.cpp)
   set_target_properties(pyulocator PROPERTIES
                         PREFIX ""
                         CXX_STANDARD 20
                         CXX_STANDARD_REQUIRED YES 
                         CXX_EXTENSIONS NO) 
   target_include_directories(pyulocator
                              PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/python/include>
                                      $<BUILD_INTERFACE:${UMPS_INCLUDE_DIR}>
                                      $<BUILD_INTERFACE:${UMPS_INCLUDE_DIR}/../>
                                      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
   target_link_libraries(pyulocator
                         PRIVATE pybind11::module pybind11::pybind11 pybind11::lto
                                 uLocator ${UMPS_LIBRARY} ${EIKONALXX_LIBRARY})
endif()

#add_executable(createSirenData src/createSirenData.cpp)
#set_target_properties(evaluateCorrections PROPERTIES
#                      CXX_STANDARD 20
#                      CXX_STANDARD_REQUIRED YES 
#                      CXX_EXTENSIONS NO) 
#target_link_libraries(createSirenData uLocator
#                      ${UMPS_LIBRARY} hdf5::hdf5 Boost::program_options
#                      GeographicLib::GeographicLib ${EIKONALXX_LIBRARY})
#target_include_directories(createSirenData
#                           PRIVATE Boost::program_options ${CMAKE_SOURCE_DIR}/include
#                                   GeographicLib::GeographicLib
#                                   ${UMPS_INCLUDE_DIR} hdf5::hdf5 ${EIKONALXX_INCLUDE_DIR})


##########################################################################################
#                                      Installation                                      #
##########################################################################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/uLocatorConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)
if (${WRAP_PYTHON})
   install(TARGETS uLocator pyulocator
           EXPORT ${PROJECT_NAME}-targets
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           COMPONENT Runtime)
else()
   install(TARGETS uLocator
           EXPORT ${PROJECT_NAME}-targets
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           COMPONENT Runtime)
endif()
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/uLocator
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "ben.baker@utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "An earthquake assocation and location utility used at UUSS.")
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
