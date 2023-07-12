add_executable(distfort main.f90)
  
add_subdirectory(distances)
add_subdirectory(linked_lists)
add_subdirectory(testing)

target_link_libraries(distfort
  PRIVATE
    distance
    linked_lists
    testing
  )