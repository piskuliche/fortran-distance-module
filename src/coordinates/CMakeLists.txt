add_library(coordinates_module "")
  
target_sources(coordinates_module
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/coordinates_module.f90
  )