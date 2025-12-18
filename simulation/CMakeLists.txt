#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PET-Simulation)
# Add option to enable/disable tests
option(ENABLE_TESTS "Enable tests" OFF)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MACRO_FILES
     "macros/*.mac"
)
file(COPY ${MACRO_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)

file(GLOB DATA_FILES
  "*.dat"
)
file(COPY ${DATA_FILES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(sim sim.cc ${sources} ${headers})
target_link_libraries(sim ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(S DEPENDS sim)

# Define some variables
set(PATH ${PROJECT_BINARY_DIR}/)
set(LIBRARY $ENV{LD_LIBRARY_PATH}/)

# Configure the config file
configure_file(
  "${PROJECT_SOURCE_DIR}/config.in"
  "${PROJECT_BINARY_DIR}/config.sh"
  @ONLY
)

if(ENABLE_TESTS)
    # Add Google Test
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})

    # Add test executable
    add_executable(test_sim test/test_sim.cc ${sources} ${headers})
    target_link_libraries(test_sim ${GTEST_BOTH_LIBRARIES} ${Geant4_LIBRARIES} pthread)

    # Add test to CTest
    add_test(NAME Geant4SimulationTest COMMAND test_sim)
else()
    message(STATUS "Tests are disabled. Use -DENABLE_TESTS=ON to enable.")
endif()
