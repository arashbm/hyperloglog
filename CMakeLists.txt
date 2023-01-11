cmake_minimum_required(VERSION 3.14)

project(hyperloglog LANGUAGES CXX VERSION 0.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(GNUInstallDirs)

add_library(murmurhash3 SHARED dep/MurmurHash3/src/MurmurHash3.cpp)
target_include_directories(
  murmurhash3
  PUBLIC
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}>/dep/MurmurHash3/include/
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(MSVC)
  target_compile_options(murmurhash3 PUBLIC /W4)
else()
  target_compile_options(murmurhash3 PUBLIC
    -Werror -Wall -Wextra -Wconversion
    -Wno-sign-conversion -Wno-implicit-fallthrough)
endif()


add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}>/include/
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(${PROJECT_NAME} INTERFACE murmurhash3)


add_executable(estimate_distribution EXCLUDE_FROM_ALL
  src/estimate_distribution.cpp)
target_link_libraries(estimate_distribution PRIVATE ${PROJECT_NAME})

add_executable(record_biases EXCLUDE_FROM_ALL
  src/record_biases.cpp)
target_link_libraries(record_biases PRIVATE ${PROJECT_NAME})

Include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.7)

FetchContent_MakeAvailable(Catch2)

add_executable(
  ${PROJECT_NAME}_tests
  EXCLUDE_FROM_ALL
  src/tests.cpp src/hll_tests.cpp)
target_link_libraries(${PROJECT_NAME}_tests PRIVATE ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_tests PRIVATE Catch2::Catch2)
if(MSVC)
  target_compile_options(${PROJECT_NAME}_tests PRIVATE /W4)
else()
  target_compile_options(${PROJECT_NAME}_tests PRIVATE
    -Werror -Wall -Wextra -Wconversion -Wsign-conversion)
endif()

install(
  TARGETS ${PROJECT_NAME} murmurhash3
  EXPORT ${PROJECT_NAME}_Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${PROJECT_NAME}-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

if(NOT INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR include/)
endif()

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
  PATH_VARS INCLUDE_INSTALL_DIR)

install(
  EXPORT ${PROJECT_NAME}_Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(
  FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
