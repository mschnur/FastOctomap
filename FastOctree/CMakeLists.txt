cmake_minimum_required(VERSION 2.8.8)
project(FastOctree C)

set(CMAKE_C_STANDARD 99)

add_executable(FastOctree main.c FastOctree.c FastOctree.h NodeStack.c NodeStack.h)