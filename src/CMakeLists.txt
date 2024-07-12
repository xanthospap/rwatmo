add_library(rwatmo)

target_sources(rwatmo
  PRIVATE
    vmf3_details.cpp
    vmf3_site_file_stream.cpp
    vmf3_site_stream.cpp
    vmf3_full_coeffs.cpp
)

target_include_directories(rwatmo
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# We need the Eigen-3 library
# target_link_libraries(geodesy PRIVATE Eigen3::Eigen)

# Install as: $PREFIX/geodesy/...
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include/
	DESTINATION rwatmo
)
