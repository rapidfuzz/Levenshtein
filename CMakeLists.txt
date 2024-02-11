cmake_minimum_required(VERSION 3.12...3.24)

cmake_policy(SET CMP0054 NEW)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(SKBUILD_LINK_LIBRARIES_KEYWORD PRIVATE)
set(Python_FIND_IMPLEMENTATIONS CPython PyPy)

set(THREADS_PREFER_PTHREAD_FLAG ON)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_OSX_DEPLOYMENT_TARGET
      "10.9"
      CACHE STRING "Minimum OS X deployment version")
endif()

project(Levenshtein LANGUAGES C CXX)

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

if(CMAKE_VERSION VERSION_LESS 3.18)
  find_package(
    Python
    COMPONENTS Interpreter Development
    REQUIRED)
else()
  set(Python_ARTIFACTS_INTERACTIVE TRUE)
  find_package(
    Python
    COMPONENTS Interpreter Development.Module
    REQUIRED)
endif()

if(CMAKE_VERSION VERSION_LESS 3.17)
  execute_process(
    COMMAND
      "${Python_EXECUTABLE}" -c
      "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX').split('.')[1])"
    OUTPUT_VARIABLE Python_SOABI
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  message(STATUS "Corrected SOABI: ${Python_SOABI}")
elseif("${Python_INTERPRETER_ID}" STREQUAL "PyPy")
  message(STATUS "PyPy SOABI: ${Python_SOABI}")
  execute_process(
    COMMAND
      "${Python_EXECUTABLE}" -c
      "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX').split('.')[1])"
    OUTPUT_VARIABLE Python_SOABI
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  message(STATUS "Corrected SOABI: ${Python_SOABI}")
endif()

find_package(rapidfuzz 3.0.0 QUIET)
if(rapidfuzz_FOUND)
  message(STATUS "Using system supplied version of rapidfuzz-cpp")
else()
  message(STATUS "Using packaged version of rapidfuzz-cpp")
  add_subdirectory(extern/rapidfuzz-cpp)
endif()

set(LEV_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/Levenshtein)

add_subdirectory(src/Levenshtein)
