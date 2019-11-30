cmake_minimum_required(VERSION 2.8.8)
project(FastOctree C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "-std=c99 -Wall -march=native -g")

add_executable(FastOctreeApp main.c FastOctree.c FastOctree.h Stack.c Stack.h)
target_link_libraries(FastOctreeApp m)

add_library(FastOctree FastOctree.c FastOctree.h Stack.c Stack.h)
target_link_libraries(FastOctree m)

install(TARGETS FastOctreeApp DESTINATION bin)
install(TARGETS FastOctree DESTINATION lib)
install(TARGETS FastOctree DESTINATION lib)
install(FILES "${CMAKE_SOURCE_DIR}/FastOctree.h" DESTINATION include)