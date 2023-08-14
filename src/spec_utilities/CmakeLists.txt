add_library(spectra_module "")
  
target_sources(spectra_module
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/spectra_module.f90
  )