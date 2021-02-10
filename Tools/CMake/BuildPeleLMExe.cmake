function(build_pelelm_exe pelelm_exe_name)

  add_executable(${pelelm_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pelelm_exe_name} PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
  endif()

  target_sources(${pelelm_exe_name}
     PRIVATE
       pelelm_prob_parm.H
       pelelm_prob.H
       pelelm_prob.cpp
  )
  
  target_include_directories(${pelelm_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/../PelePhysics)
  set(PELE_PHYSICS_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/PelePhysics/${pelelm_exe_name})

  set(IAMR_SRC_DIR ${CMAKE_SOURCE_DIR}/../IAMR/Source)
  set(IAMR_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/IAMR/Source/${pelelm_exe_name})

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${pelelm_exe_name})

  include(SetPeleLMCompileFlags)

  set(PELELM_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport/${PELELM_TRANSPORT_MODEL}")
  target_sources(${pelelm_exe_name} PRIVATE
                 ${PELELM_TRANSPORT_DIR}/Transport.H
                 ${PELELM_TRANSPORT_DIR}/Transport.cpp
                 ${PELELM_TRANSPORT_DIR}/TransportParams.H)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELELM_TRANSPORT_DIR})

  set(PELELM_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos/${PELELM_EOS_MODEL}")
  target_sources(${pelelm_exe_name} PRIVATE
                 ${PELELM_EOS_DIR}/EOS.cpp
                 ${PELELM_EOS_DIR}/EOS.H)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELELM_EOS_DIR})

  set(PELELM_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Mechanism/Models/${PELELM_CHEMISTRY_MODEL}")
  target_sources(${pelelm_exe_name} PRIVATE
                 ${PELELM_MECHANISM_DIR}/chemistry_file.H
                 ${PELELM_MECHANISM_DIR}/mechanism.cpp
                 ${PELELM_MECHANISM_DIR}/mechanism.h)
  # Avoid warnings from mechanism.cpp for now
  if(NOT PELELM_ENABLE_CUDA)
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
        list(APPEND MY_CXX_FLAGS "-Wno-unreachable-code"
                                 "-Wno-null-dereference"
                                 "-Wno-float-conversion"
                                 "-Wno-shadow"
                                 "-Wno-overloaded-virtual")
      endif()
      if(CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|AppleClang)$")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-variable")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-parameter")
        list(APPEND MY_CXX_FLAGS "-Wno-vla-extension")
        list(APPEND MY_CXX_FLAGS "-Wno-zero-length-array")
      elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-variable")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-parameter")
        list(APPEND MY_CXX_FLAGS "-Wno-vla")
        list(APPEND MY_CXX_FLAGS "-Wno-pedantic")
      endif()
    endif()
  endif()
  separate_arguments(MY_CXX_FLAGS)
  set_source_files_properties(${PELELM_MECHANISM_DIR}/mechanism.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  set_source_files_properties(${PELELM_MECHANISM_DIR}/chemistry_file.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  set_source_files_properties(${PELELM_MECHANISM_DIR}/mechanism.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELELM_MECHANISM_DIR})
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Evaluation)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Utility)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Mechanism/Models)
  
  if(PELELM_ENABLE_SUNDIALS)
    if(PELELM_ENABLE_CUDA)
      set(DEVICE GPU)
    else()
      set(DEVICE CPU)
    endif()
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_SUNDIALS_PP)
    target_sources(${pelelm_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/cvode/reactor.cpp
                                             ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/cvode/reactor.h)
    set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/cvode/reactor.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
    set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/cvode/reactor.h PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
    target_link_libraries(${pelelm_exe_name} PRIVATE sundials_cvode)
    target_include_directories(${pelelm_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE})
    target_include_directories(${pelelm_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/cvode)
    if(PELELM_ENABLE_CUDA)
      target_sources(${pelelm_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/AMReX_SUNMemory.cpp
                                                ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/AMReX_SUNMemory.H)
      target_link_libraries(${pelelm_exe_name} PRIVATE sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
    endif()
  endif()
  
  target_sources(${pelelm_exe_name}
     PRIVATE
       ${SRC_DIR}/ArrayViewEXT.H
       ${SRC_DIR}/BoxLib_Data_Dump.H
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/PeleLM.H
       ${SRC_DIR}/PeleLM.cpp
       ${SRC_DIR}/PeleLM_K.H
       ${SRC_DIR}/PeleLM_bcfill.cpp
       ${SRC_DIR}/PeleLM_derive.H
       ${SRC_DIR}/PeleLM_derive.cpp
       ${SRC_DIR}/PeleLM_parm.H
       ${SRC_DIR}/PeleLM_setup.cpp
       ${SRC_DIR}/Prob_F.H
  )

  target_sources(${pelelm_exe_name}
     PRIVATE
       ${IAMR_SRC_DIR}/SyncRegister.cpp
       ${IAMR_SRC_DIR}/NS_init_eb2.cpp
       ${IAMR_SRC_DIR}/NS_getForce.cpp
       ${IAMR_SRC_DIR}/OutFlowBC.cpp
       ${IAMR_SRC_DIR}/FluxBoxes.cpp
       ${IAMR_SRC_DIR}/OutFlowBC.H
       ${IAMR_SRC_DIR}/SyncRegister.H
       ${IAMR_SRC_DIR}/RegType.H
       ${IAMR_SRC_DIR}/MacOpMacDrivers.H
       ${IAMR_SRC_DIR}/iamr_constants.H
       ${IAMR_SRC_DIR}/NavierStokesBase.cpp
       ${IAMR_SRC_DIR}/Projection.cpp
       ${IAMR_SRC_DIR}/MacProj.cpp
       ${IAMR_SRC_DIR}/Diffusion.cpp
       ${IAMR_SRC_DIR}/NS_LES.cpp
       ${IAMR_SRC_DIR}/NS_derive.cpp
       ${IAMR_SRC_DIR}/NS_average.cpp
       ${IAMR_SRC_DIR}/NS_derive.H
       ${IAMR_SRC_DIR}/Projection.H
       ${IAMR_SRC_DIR}/MacProj.H
       ${IAMR_SRC_DIR}/Diffusion.H
       ${IAMR_SRC_DIR}/NavierStokesBase.H
       ${IAMR_SRC_DIR}/FluxBoxes.H
       ${IAMR_SRC_DIR}/DERIVE_F.H
       ${IAMR_SRC_DIR}/PROJECTION_F.H
       ${IAMR_SRC_DIR}/NAVIERSTOKES_F.H
       ${IAMR_SRC_DIR}/NS_error_F.H
       ${IAMR_SRC_DIR}/MLMG_Mac.cpp
       ${IAMR_SRC_DIR}/NS_util.cpp
       ${IAMR_SRC_DIR}/MACPROJ_F.H
       ${IAMR_SRC_DIR}/NS_util.H

       ${IAMR_SRC_DIR}/Godunov/iamr_godunov.H
       ${IAMR_SRC_DIR}/Godunov/iamr_godunov_K.H
       ${IAMR_SRC_DIR}/Godunov/iamr_constants_godunov.H
       ${IAMR_SRC_DIR}/Godunov/iamr_godunov.cpp
       ${IAMR_SRC_DIR}/Godunov/iamr_extrap_vel_to_faces_godunov_${PELELM_DIM}D.cpp
       ${IAMR_SRC_DIR}/Godunov/iamr_edge_state_godunov_${PELELM_DIM}D.cpp
       ${IAMR_SRC_DIR}/Godunov/iamr_slopes_godunov_K.H
       ${IAMR_SRC_DIR}/Godunov/iamr_plm_godunov.H
       ${IAMR_SRC_DIR}/Godunov/iamr_plm_godunov.cpp
       ${IAMR_SRC_DIR}/Godunov/iamr_ppm_godunov.H
       ${IAMR_SRC_DIR}/Godunov/iamr_ppm_godunov.cpp

       # Dimension-specific files
       ${IAMR_SRC_DIR}/Src_${PELELM_DIM}d/NAVIERSTOKES_${PELELM_DIM}D.F90
       ${IAMR_SRC_DIR}/Src_${PELELM_DIM}d/PROJECTION_${PELELM_DIM}D.F90
       ${IAMR_SRC_DIR}/Src_${PELELM_DIM}d/MACPROJ_${PELELM_DIM}D.F90
  )

  if(PELELM_ENABLE_AMREX_EB)
    target_sources(${pelelm_exe_name}
       PRIVATE
         ${IAMR_SRC_DIR}/MOL/iamr_mol.H
         ${IAMR_SRC_DIR}/MOL/iamr_mol.cpp
         ${IAMR_SRC_DIR}/MOL/iamr_edge_state_mol.cpp
         ${IAMR_SRC_DIR}/MOL/iamr_extrap_vel_to_faces_mol.cpp
         ${IAMR_SRC_DIR}/MOL/iamr_edge_state_mol_K.H
         ${IAMR_SRC_DIR}/MOL/iamr_predict_vel_on_faces_mol.cpp
         ${IAMR_SRC_DIR}/MOL/iamr_eb_predict_vel_on_faces_mol.cpp
         ${IAMR_SRC_DIR}/MOL/iamr_eb_edge_state_mol_K.H
         #${IAMR_SRC_DIR}/MOL/iamr_eb_slopes_mol_K.H
         #${IAMR_SRC_DIR}/MOL/iamr_slopes_K.H
    )
  endif()

  target_sources(${pelelm_exe_name}
     PRIVATE
       ${PELE_PHYSICS_SRC_DIR}/Utility/pmf.cpp
       ${PELE_PHYSICS_SRC_DIR}/Utility/pmf_data.cpp
       ${PELE_PHYSICS_SRC_DIR}/Utility/pmf.H
       ${PELE_PHYSICS_SRC_DIR}/Utility/pmf_data.H
  )

  if(NOT "${pelelm_exe_name}" STREQUAL "PeleLM-UnitTests")
    target_sources(${pelelm_exe_name}
       PRIVATE
         ${IAMR_SRC_DIR}/main.cpp
    )
  endif()
  
  include(AMReXBuildInfo)
  generate_buildinfo(${pelelm_exe_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(PELELM_ENABLE_MASA)
    target_compile_definitions(${pelelm_exe_name} PRIVATE PELELM_USE_MASA)
    #target_sources(${pelelm_exe_name} PRIVATE ${SRC_DIR}/MMS.cpp)
    target_link_libraries(${pelelm_exe_name} PRIVATE MASA::MASA)
    #if(PELELM_ENABLE_FPE_TRAP_FOR_TESTS)
    #  set_source_files_properties(${SRC_DIR}/PeleLM.cpp PROPERTIES COMPILE_DEFINITIONS PELELM_ENABLE_FPE_TRAP)
    #endif()
  endif()

  if(PELELM_ENABLE_MPI)
    target_link_libraries(${pelelm_exe_name} PRIVATE $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #PeleLM include directories
  target_include_directories(${pelelm_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pelelm_exe_name} PRIVATE ${CMAKE_BINARY_DIR})

  #IAMR include directories
  target_include_directories(${pelelm_exe_name} PRIVATE ${IAMR_SRC_DIR})
  target_include_directories(${pelelm_exe_name} PRIVATE ${IAMR_SRC_DIR}/Godunov)
  if(PELELM_ENABLE_AMREX_EB)
    target_include_directories(${pelelm_exe_name} PRIVATE ${IAMR_SRC_DIR}/MOL)
  endif()

  #Keep our Fortran module files confined to a unique directory for each executable 
  set_target_properties(${pelelm_exe_name} PROPERTIES Fortran_MODULE_DIRECTORY
                       "${CMAKE_BINARY_DIR}/fortran_modules/${pelelm_exe_name}_fortran_modules")
  target_include_directories(${pelelm_exe_name} PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules/${pelelm_exe_name}_fortran_modules)

  #Link to amrex library
  target_link_libraries(${pelelm_exe_name} PRIVATE AMReX::amrex)

  if(PELELM_ENABLE_CUDA)
    set(plmtargets "${pelelm_exe_name}")
    foreach(tgt IN LISTS plmtargets)
      get_target_property(PELELM_SOURCES ${tgt} SOURCES)
      list(FILTER PELELM_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${PELELM_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${pelelm_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${pelelm_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()
 
  #Define what we want to be installed during a make install 
  install(TARGETS ${pelelm_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
