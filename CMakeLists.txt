cmake_minimum_required(VERSION 3.1)

project(
  rwatmo
  VERSION 1.0.1
  DESCRIPTION "Radio-Wave Atmosphere"
  LANGUAGES CXX)


# Ensure required libraries are available
find_package(Eigen3   REQUIRED)
find_package(iers     REQUIRED)
find_package(geodesy  REQUIRED)
find_package(datetime REQUIRED)

# Pass the library dependencies to subdirectories
set(PROJECT_DEPENDENCIES Eigen3::Eigen iers geodesy datetime)

# Define an option for building tests (defaults to ON)
option(BUILD_TESTING "Enable building of tests" ON)

# compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
add_compile_options(-Wall
  -Wextra
  -Werror
  -pedantic
  -W
  -Wshadow
  $<$<CONFIG:Release>:-O2>
  $<$<CONFIG:Release>:-march=native>
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Debug>:-pg>
  $<$<CONFIG:Debug>:-Wdisabled-optimization>
  )
add_compile_definitions(
  $<$<CONFIG:Debug>:DEBUG>
)

# the library and includes
add_library(rwatmo)
target_include_directories(rwatmo
PUBLIC
  # Public headers during build time
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  # Public headers during install time
  $<INSTALL_INTERFACE:include/rwatmo>
)

# The library
add_subdirectory(src/vmf)
add_subdirectory(src/nrlmsise00)
add_subdirectory(src/math)
add_subdirectory(src/spaceweather)

# The tests
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(test/)
  enable_testing()
endif()

# install library
install(TARGETS rwatmo
        EXPORT rwatmoTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)

# Install public headers at: $PREFIX/rwatmo/...
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
	DESTINATION include/rwatmo
)

install(EXPORT rwatmoTargets
        FILE rwatmoTargets.cmake
        NAMESPACE dso::
        DESTINATION lib/cmake
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "rwatmoConfigVersion.cmake"
  VERSION ${rwatmo_version}
  COMPATIBILITY AnyNewerVersion
)

install(FILES "rwatmoConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/rwatmoConfigVersion.cmake"
  DESTINATION lib/cmake
)
