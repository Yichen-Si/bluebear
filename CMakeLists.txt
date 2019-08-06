cmake_minimum_required (VERSION 2.6)
project(ibs0phase)

include_directories(${CMAKE_SOURCE_DIR}/../htslib/include)

execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)
set(VERSION "1.0.0")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 -DVERSION='\"${VERSION}\"' -DUSER='\"${USER}\"' -DDATE='\"${DATE}\"' -D__ZLIB_AVAILABLE__")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(src)
file(GLOB_RECURSE ibs0phase_SOURCES "src/*.cpp" "src/*.c")
file(GLOB_RECURSE ibs0phase_HEADERS "src/*.h")

add_executable(bluebear ${ibs0phase_SOURCES})
target_link_libraries(bluebear ${CMAKE_SOURCE_DIR}/../htslib/libhts.a z dl pthread lzma z bz2 curl crypto)

