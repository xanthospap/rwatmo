target_sources(rwatmo
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/vmf/load_sites_for_epoch.cpp
    ${CMAKE_SOURCE_DIR}/src/vmf/load_vmf3gr_grid_data.cpp
    ${CMAKE_SOURCE_DIR}/src/vmf/vmf3_static_storage.cpp
    ${CMAKE_SOURCE_DIR}/src/vmf/vmf3_compute_cell.cpp
)
