target_sources(rwatmo
  PRIVATE
    #${CMAKE_SOURCE_DIR}/src/vmf/vmf3_details.cpp
    #${CMAKE_SOURCE_DIR}/src/vmf/vmf3_full_coeffs.cpp
    #${CMAKE_SOURCE_DIR}/src/vmf/vmf3_site_file_stream.cpp
    #${CMAKE_SOURCE_DIR}/src/vmf/vmf3_site_stream.cpp
    ${CMAKE_SOURCE_DIR}/src/vmf/load_gr_grid_data.cpp
)
