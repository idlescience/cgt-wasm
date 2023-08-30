set(
        MODULE_FILES
        src/cgt_module.cpp
)

add_executable (cgt_module ${MODULE_FILES})
target_link_libraries(cgt_module cgt embind)

set(MODULE_INCLUDE_DIR include)

target_include_directories(cgt_module PUBLIC ${MODULE_INCLUDE_DIR})

set_target_properties(
  cgt_module
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LINK_FLAGS "-s ENVIRONMENT=web -s SINGLE_FILE=1 -s MODULARIZE=1 -s EXPORT_ES6=1 -s EXPORT_NAME=cgt_module --bind"
)
