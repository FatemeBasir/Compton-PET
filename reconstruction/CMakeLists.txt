cmake_minimum_required(VERSION 3.10)

# Project name
project(PET-Digitization)

# Find the ROOT package
find_package(ROOT REQUIRED)

# Include ROOT directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Create an executable
add_executable(reco reco.cc ${sources} ${headers})

# Link the ROOT libraries
target_link_libraries(reco ${ROOT_LIBRARIES})
