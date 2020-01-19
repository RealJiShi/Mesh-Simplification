cmake_minimum_required(VERSION 3.4.1)

# settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-rtti -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# include
set(CMAKE_INCLUDE_CURRENT_DIR)

# source
aux_source_directory(./ SRC)

#
add_executable(MeshSimp ${SRC})
