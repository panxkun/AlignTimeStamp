cmake_minimum_required(VERSION 3.0)

project(align_timestamp)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++14 -O3" )

# 寻找Ceres库并添加它的头文件
find_package( Ceres REQUIRED )
include_directories( ${CERES_INCLUDE_DIRS} )

include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/)

# add_executable(leastsquare leastsquare.cpp)
# target_link_libraries(leastsquare ${CERES_LIBRARIES})

add_executable(ncc ncc.cpp)
target_link_libraries(ncc ${CERES_LIBRARIES})