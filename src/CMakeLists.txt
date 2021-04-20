# Make an automatic library - will be static or dynamic based on user setting
add_library(shapelet_transform shapelet_transform.cpp)

target_include_directories(shapelet_transform PUBLIC ../include)

# This depends on OpenCL 
include(FindOpenCL)

target_link_libraries(shapelet_transform PRIVATE OpenCL::OpenCL)

# All users of this library will need at least C++11
target_compile_features(shapelet_transform PUBLIC cxx_std_11)