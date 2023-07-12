add_library(distance_module "")
  
target_sources(distance_module
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/distance_module.f90
  )