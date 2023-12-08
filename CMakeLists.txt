cmake_minimum_required(VERSION 3.12)

project(2D_Kuramoto-Sivashinsky)
# set compiler
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)


# set variables
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SOURCE_FILES 
  ${SOURCE_DIR}/main.cpp
  ${SOURCE_DIR}/fft.cpp
  ${SOURCE_DIR}/setup.cpp
  ${SOURCE_DIR}/vector.cpp
)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(TARGET_FILE prova_plot)

# find packages
find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)
find_package(matplotlibcpp17 REQUIRED)

add_executable(${TARGET_FILE} ${SOURCE_FILES})
target_include_directories(${TARGET_FILE} PUBLIC ${Python3_INCLUDE_DIRS})
target_include_directories(${TARGET_FILE} PRIVATE ${INCLUDE_DIR})
target_link_libraries(${TARGET_FILE} ${Python3_LIBRARIES} pybind11::embed matplotlibcpp17::matplotlibcpp17)