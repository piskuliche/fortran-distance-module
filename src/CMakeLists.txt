add_executable(distfort main.f90)
add_executable(efield_fort efield_main.f90)

add_subdirectory(coordinates)
add_subdirectory(distances)
add_subdirectory(linked_lists)
add_subdirectory(testing)
add_subdirectory(efield_utilities)

target_link_libraries(coordinates_module)

target_link_libraries(distance_module
  PUBLIC
    MPI
  )

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
    coordinates_module
  PUBLIC
    MPI
  )

target_link_libraries(distfort
  PRIVATE
    testing
  PUBLIC
    MPI
  )

target_link_libraries(efield_fort
  PRIVATE
    distance_module
    coordinates_module
    efield_module
  PUBLIC
    MPI
    gmxfort
  )
