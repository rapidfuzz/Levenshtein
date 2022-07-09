
add_cython_target(levenshtein_cpp CXX)
add_library(levenshtein_cpp MODULE ${levenshtein_cpp} ${LEV_BASE_DIR}/Levenshtein-c/_levenshtein.cpp)
target_compile_features(levenshtein_cpp PUBLIC cxx_std_14)
target_include_directories(levenshtein_cpp PRIVATE ${LEV_BASE_DIR}/Levenshtein-c)
target_link_libraries(levenshtein_cpp PRIVATE rapidfuzz::rapidfuzz)
python_extension_module(levenshtein_cpp)
install(TARGETS levenshtein_cpp LIBRARY DESTINATION src/Levenshtein)
