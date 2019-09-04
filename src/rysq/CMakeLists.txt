add_library(librysq SHARED rysq.cc)

include_directories(${PROJECT_SOURCE_DIR}/src)

include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${PROJECT_SOURCE_DIR}/lib)

if (CUDA_FOUND)
  add_subdirectory(cuda)
endif()

if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  add_executable(validate-executable  EXCLUDE_FROM_ALL validate.cc)
  target_link_libraries(validate-executable librysq)
  add_custom_target(validate
    DEPENDS validate-executable
    COMMAND validate-executable
  )
endif()


add_executable(timings-executable  EXCLUDE_FROM_ALL timings.cc)
target_link_libraries(timings-executable librysq)

add_custom_target(timings
  DEPENDS timings-executable
  COMMAND timings-executable
  )
