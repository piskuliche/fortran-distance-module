add_library(efield_module "")
  
target_sources(efield_module
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/efield_module.f90
  )