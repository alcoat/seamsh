cmake_minimum_required (VERSION 2.8.11)
project(seamesh C CXX)
add_library(seamesh SHARED seamesh.c shp.cpp)
target_link_libraries(seamesh PUBLIC proj)
