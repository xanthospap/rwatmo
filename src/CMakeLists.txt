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
# target_link_libraries(rwatmo PRIVATE Eigen3::Eigen)
target_link_libraries(rwatmo PRIVATE geodesy)
target_link_libraries(rwatmo PRIVATE datetime)

# Install headers at: $PREFIX/rwatmo/...
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include/
	DESTINATION include/rwatmo
)

# install library
install(TARGETS rwatmo 
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin)
