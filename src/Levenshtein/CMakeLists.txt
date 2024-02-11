function(create_cython_target _name)
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/${_name}.cxx)
    set(${_name}
        ${CMAKE_CURRENT_LIST_DIR}/${_name}.cxx
        PARENT_SCOPE)
  else()
    find_package(Cython REQUIRED)
    add_cython_target(${_name} CXX)
    set(${_name}
        ${_name}
        PARENT_SCOPE)
  endif()
endfunction(create_cython_target)

function(rf_add_library name)
  if(CMAKE_VERSION VERSION_LESS 3.17)
    python_add_library(${name} MODULE ${ARGV})
    get_property(
      suffix
      TARGET ${name}
      PROPERTY SUFFIX)
    if(NOT suffix)
      set(suffix "${CMAKE_SHARED_MODULE_SUFFIX}")
    endif()
    set_property(TARGET ${name} PROPERTY SUFFIX ".${Python_SOABI}${suffix}")
  else()
    python_add_library(${name} MODULE WITH_SOABI ${ARGV})
  endif()
endfunction(rf_add_library)

create_cython_target(levenshtein_cpp)
rf_add_library(levenshtein_cpp ${levenshtein_cpp}
               ${LEV_BASE_DIR}/Levenshtein-c/_levenshtein.cpp)
target_compile_features(levenshtein_cpp PUBLIC cxx_std_17)
target_include_directories(levenshtein_cpp
                           PRIVATE ${LEV_BASE_DIR}/Levenshtein-c)
target_link_libraries(levenshtein_cpp PRIVATE rapidfuzz::rapidfuzz)
install(TARGETS levenshtein_cpp LIBRARY DESTINATION src/Levenshtein)
