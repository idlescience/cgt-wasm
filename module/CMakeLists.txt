set(
        MODULE_FILES
        src/nucleolus_module.cpp
)

add_executable (nucleolus_module ${MODULE_FILES})
target_link_libraries(nucleolus_module nucleolus embind)

set(MODULE_INCLUDE_DIR include)

target_include_directories(nucleolus_module PUBLIC ${MODULE_INCLUDE_DIR})

set_target_properties(
  nucleolus_module
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LINK_FLAGS "-s ENVIRONMENT=web -s SINGLE_FILE=1 -s MODULARIZE=1 -s EXPORT_ES6=1 -s EXPORT_NAME=nucleolus_module --bind"
)
