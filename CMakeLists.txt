cmake_minimum_required(VERSION 2.6)
project(zshrink)

set(CMAKE_C_FLAGS "-g -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")
add_definitions(-DDEBUG)
add_subdirectory(zshrink)
