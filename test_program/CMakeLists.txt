cmake_minimum_required(VERSION 2.8.8)
project(test_program)

if(NOT DEFINED OCTOMAP_DEPLOY_DIR)
	message(FATAL_ERROR "OCTOMAP_DEPLOY_DIR must be defined")
endif()

if(NOT DEFINED TEST_EXECUTABLE_NAME)
	message(FATAL_ERROR "TEST_EXECUTABLE_NAME must be defined")
endif()

if(NOT IS_ABSOLUTE ${OCTOMAP_DEPLOY_DIR})
    set(OCTOMAP_DEPLOY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${OCTOMAP_DEPLOY_DIR}")
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

find_package(octomap REQUIRED PATHS ${OCTOMAP_DEPLOY_DIR} NO_DEFAULT_PATH)

add_executable(${TEST_EXECUTABLE_NAME} test_main.cpp)

target_include_directories(${TEST_EXECUTABLE_NAME} PUBLIC ${OCTOMAP_INCLUDE_DIRS})
target_link_libraries(${TEST_EXECUTABLE_NAME} 
	PUBLIC ${OCTOMAP_DEPLOY_DIR}/lib/liboctomap.a
    PUBLIC ${OCTOMAP_DEPLOY_DIR}/lib/liboctomath.a
	)

set_target_properties(${TEST_EXECUTABLE_NAME} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -g")
#set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -g -save-temps -fverbose-asm -S")
#set(CMAKE_CXX_FLAGS "-std=c++11 -DUSE_REVELLES_RAY_TRACE_MOD_NODE=1")
add_compile_options(-std=c++11)



install(TARGETS ${TEST_EXECUTABLE_NAME} DESTINATION bin)