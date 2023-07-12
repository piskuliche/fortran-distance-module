add_library(testing "")
  
target_sources(testing
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/testing_module.f90
  )