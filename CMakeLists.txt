cmake_minimum_required(VERSION 2.6)
project(zshrink)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

set(CMAKE_C_FLAGS "-Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")
set(CMAKE_C_FLAGS_DEBUG "-DDEBUG")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_RELEASE "-DRELEASE")
add_subdirectory(zshrink)
