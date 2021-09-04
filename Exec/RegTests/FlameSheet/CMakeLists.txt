#set(PELELM_ENABLE_PARTICLES OFF)
set(PELELM_EOS_MODEL Fuego)
set(PELELM_CHEMISTRY_MODEL drm19)
set(PELELM_TRANSPORT_MODEL Simple)
set(PELELM_REACTOR_MODEL cvode)

#Generic stuff for every executable
get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(pelelm_exe_name PeleLM-${DIR_NAME})
include(BuildPeleLMExe)
build_pelelm_exe(${pelelm_exe_name})
