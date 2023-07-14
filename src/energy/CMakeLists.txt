add_library(energy "")
  
target_sources(energy
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/energy_module.f90
  )