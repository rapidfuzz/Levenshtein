cmake_minimum_required(VERSION 3.12.0)

set(CMAKE_BUILD_TYPE "RelWithDebInfo")

cmake_policy(SET CMP0054 NEW)
set(SKBUILD_LINK_LIBRARIES_KEYWORD PRIVATE)

set(THREADS_PREFER_PTHREAD_FLAG ON)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")
endif()

project(Levenshtein LANGUAGES C CXX)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

find_package(PythonExtensions REQUIRED)
find_package(Python COMPONENTS Interpreter Development)

find_package(rapidfuzz 1.0.4 QUIET)
if (rapidfuzz_FOUND)
    message("Using system supplied version of rapidfuzz-cpp")
else()
    message("Using packaged version of rapidfuzz-cpp")
    add_subdirectory(extern/rapidfuzz-cpp)
endif()

set(LEV_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/Levenshtein)

add_subdirectory(src/Levenshtein)
