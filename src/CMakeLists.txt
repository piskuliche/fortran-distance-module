add_executable(distfort main.f90)


add_subdirectory(distances)
add_subdirectory(linked_lists)
add_subdirectory(testing)

target_link_libraries(distance_module)

target_link_libraries(linked_lists
  PRIVATE
    distance_module
  PUBLIC
    MPI
  )

target_link_libraries(testing 
  PRIVATE
    distance_module
    linked_lists
  PUBLIC
    MPI
  )

target_link_libraries(distfort
  PRIVATE
    testing
  PUBLIC
    MPI
  )
