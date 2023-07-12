add_library(distance "")
  
target_sources(distance
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/distance_module.f90
  )