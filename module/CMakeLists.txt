set(
        MODULE_FILES
        src/bnf_module.cpp
)

add_executable (module ${MODULE_FILES})
target_link_libraries(module nucleolus embind)

set(MODULE_INCLUDE_DIR include)

target_include_directories(module PUBLIC ${MODULE_INCLUDE_DIR})

set_target_properties(
  module
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
)