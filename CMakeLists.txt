cmake_minimum_required(VERSION 2.8)


set(CMAKE_BUILD_TYPE Debug CACHE STRING "Debug | Release")
set(DEBUG ON CACHE BOOL "for the internal usage - not the compiler debug options")

if(CMAKE_BUILD_TYPE MATCHES Release)
  set(DEBUG OFF CACHE BOOL "for the internal usage - not the compiler debug options" FORCE)
endif(CMAKE_BUILD_TYPE MATCHES Release)


# --- define and connect some important directories ---
set(INSTALL_DIR $ENV{HOME}/install)
set(INCLUDE_DIR ${INSTALL_DIR}/include)
set(LIB_DIR     ${INSTALL_DIR}/lib)
include_directories(${INCLUDE_DIR})
link_directories(${LIB_DIR})
# -----------------------------------------------------


# --- project ---
#set(CMAKE_CXX_COMPILER mpicxx)
project(fem)
add_definitions(-std=c++0x)
# ---------------


# --- boost ---
#find_library(BOOST_PO_LIB libboost_program_options.so ${LIB_DIR})
#if(NOT BOOST_PO_LIB)
#  message(FATAL_ERROR "boost program_options library wasn't found!")
#endif(NOT BOOST_PO_LIB)
#find_library(BOOST_TIMER_LIB libboost_timer.so ${LIB_DIR})
#if(NOT BOOST_TIMER_LIB)
#  message(FATAL_ERROR "boost timer library wasn't found!")
#endif(NOT BOOST_TIMER_LIB)
set(Boost_INCLUDE_DIR ${INCLUDE_DIR})
set(Boost_LIBRARY_DIR ${LIB_DIR})
find_package(Boost 1.55 COMPONENTS program_options timer system filesystem REQUIRED)
# -------------


# --- google test ---
#find_library(GTEST_LIB gtest ${LIB_DIR})
#if(NOT GTEST_LIB)
#  message(FATAL_ERROR "GTest library was not found!")
#endif(NOT GTEST_LIB)
# -------------------


# --- petsc ---
find_library(PETSC_LIB petsc ${LIB_DIR})
if(NOT PETSC_LIB)
  message(FATAL_ERROR "PETSc library was not found!")
endif(NOT PETSC_LIB)
# -------------


# --- mpi ---
#set(MPI_LIB libmpi.so ${LIB_DIR})
#if(NOT MPI_LIB)
#  message(FATAL_ERROR "mpi library was not found!")
#endif(NOT MPI_LIB)
# -----------


# --- 64-bitness ---
set(HAVE_64BIT_SIZE_T OFF CACHE INTERNAL "")
include(CheckTypeSize)
check_type_size("void*" SIZEOF_VOID_P)
if(SIZEOF_VOID_P EQUAL 8)
  set(HAVE_64BIT_SIZE_T ON CACHE INTERNAL "" FORCE)
endif(SIZEOF_VOID_P EQUAL 8)
# ------------------


# --- config file ---
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_SOURCE_DIR}/headers/config.h")
# -------------------


# --- headers and sources ---
aux_source_directory(${PROJECT_SOURCE_DIR}/sources SRC_LIST) # all .cpp files
include_directories(${PROJECT_SOURCE_DIR}/headers)
FILE(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/headers/*.h")      # all .h files
# ---------------------------


# ------
message("build type   " ${CMAKE_BUILD_TYPE})
#message("gtest lib    " ${GTEST_LIB})
message("boost libs   " ${Boost_LIBRARIES})
message("petsc lib    " ${PETSC_LIB})
#message("mpi lib      " ${MPI_LIB})
message("SRC_LIST     " ${SRC_LIST})
message("HDR_LIST     " ${HDR_LIST})
# ------


#if(CMAKE_BUILD_TYPE MATCHES Debug)
#  set(LIB_NAME "femd")
#else(CMAKE_BUILD_TYPE MATCHES Debug)
#  set(LIB_NAME "fem")
#endif(CMAKE_BUILD_TYPE MATCHES Debug)


add_library(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${PETSC_LIB}) # ${GTEST_LIB} ${PETSC_LIB} ${MPI_LIB})


install(TARGETS ${PROJECT_NAME} DESTINATION ${LIB_DIR})
install(FILES ${HDR_LIST} DESTINATION ${INCLUDE_DIR}/fem)
