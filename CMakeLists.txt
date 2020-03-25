cmake_minimum_required (VERSION 2.8.11)
project(seamesh C)
add_library(seamesh SHARED seamesh.c)
target_link_libraries(seamesh PUBLIC proj)
