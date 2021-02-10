#set(PELELM_ENABLE_REACTIONS OFF)
set(PELELM_ENABLE_PARTICLES OFF)
set(PELELM_EOS_MODEL Fuego)
set(PELELM_CHEMISTRY_MODEL drm19)
set(PELELM_TRANSPORT_MODEL Simple)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(pelelm_exe_name PeleLM-${DIR_NAME})
include(BuildPeleLMExe)
build_pelelm_exe(${pelelm_exe_name})
target_sources(${pelelm_exe_name}
  PUBLIC
  unit-tests-main.cpp
  test-config.cpp
  )

if(PELELM_ENABLE_CUDA)
  set_source_files_properties(unit-tests-main.cpp test-config.cpp PROPERTIES LANGUAGE CUDA)
endif()

target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/../googletest/googletest/include)
target_link_libraries(${pelelm_exe_name} PRIVATE gtest)
