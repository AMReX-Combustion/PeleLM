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

  set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
  set(PELE_PHYSICS_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/PelePhysics/${pelelm_exe_name})

  set(IAMR_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/IAMR/Source)
  set(IAMR_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/IAMR/Source/${pelelm_exe_name})

  set(AMREX_HYDRO_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/AMReX-Hydro/)
  set(AMREX_HYDRO_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/AMReX-Hydro/${pelelm_exe_name})

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${pelelm_exe_name})

  include(SetPeleLMCompileFlags)

  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE "${PELE_PHYSICS_SRC_DIR}/Source")

  set(PELELM_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport")
  target_sources(${pelelm_exe_name} PRIVATE
                 ${PELELM_TRANSPORT_DIR}/Transport.H
                 ${PELELM_TRANSPORT_DIR}/Transport.cpp
                 ${PELELM_TRANSPORT_DIR}/TransportParams.H
                 ${PELELM_TRANSPORT_DIR}/TransportTypes.H
                 ${PELELM_TRANSPORT_DIR}/Simple.H
                 ${PELELM_TRANSPORT_DIR}/Sutherland.H
                 ${PELELM_TRANSPORT_DIR}/Constant.H)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELELM_TRANSPORT_DIR})
  if (PELELM_TRANSPORT_MODEL STREQUAL "Simple")
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_SIMPLE_TRANSPORT)
  elseif(PELELM_TRANSPORT_MODEL STREQUAL "EGLib")
    target_compile_definitions(${pelelm_exe_name} PRIVATE EGLIB_TRANSPORT)
  elseif(PELELM_TRANSPORT_MODEL STREQUAL "Constant")
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_CONSTANT_TRANSPORT)
  elseif(PELELM_TRANSPORT_MODEL STREQUAL "Sutherland")
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_SUTHERLAND_TRANSPORT)
  endif()

  set(PELELM_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos")
  target_sources(${pelelm_exe_name} PRIVATE
                 ${PELELM_EOS_DIR}/EOS.cpp
                 ${PELELM_EOS_DIR}/EOS.H
                 ${PELELM_EOS_DIR}/Fuego.H
                 ${PELELM_EOS_DIR}/GammaLaw.H
                 ${PELELM_EOS_DIR}/SRK.H)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELELM_EOS_DIR})
  if (PELELM_EOS_MODEL STREQUAL "Fuego")
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_FUEGO_EOS)
  elseif(PELELM_EOS_MODEL STREQUAL "GammaLaw")
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_GAMMALAW_EOS)
  elseif(PELELM_EOS_MODEL STREQUAL "Soave-Redlich-Kwong")
    target_compile_definitions(${pelelm_exe_name} PRIVATE USE_SRK_EOS)
  endif()

  set(PELELM_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Mechanism/Models/${PELELM_CHEMISTRY_MODEL}")
  target_sources(${pelelm_exe_name} PRIVATE
                 ${PELELM_MECHANISM_DIR}/mechanism.cpp
                 ${PELELM_MECHANISM_DIR}/mechanism.H)
  # Set PeleMP flags
  set(PELEMP_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PeleMP/Source)
  if(PELELM_ENABLE_PARTICLES)
    target_sources(${pelelm_exe_name}
      PRIVATE
	SprayParticlesInitInsert.cpp
    )
    target_compile_definitions(${pelelm_exe_name} PRIVATE SPRAY_FUEL_NUM=${PELEMP_SPRAY_FUEL_NUM})
    target_compile_definitions(${pelelm_exe_name} PRIVATE SPRAY_PELE_LM)
    target_sources(${pelelm_exe_name} PRIVATE
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayParticles.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayParticles.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayFuelData.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayInterpolation.H
                   ${PELEMP_SRC_DIR}/PP_Spray/Drag.H
                   ${PELEMP_SRC_DIR}/PP_Spray/WallFunctions.H)
    target_include_directories(${pelelm_exe_name} PRIVATE ${PELEMP_SRC_DIR}/PP_Spray)
  endif()
  if(PELELM_ENABLE_SOOT AND PELELM_SOOT_MODEL)
    target_compile_definitions(${pelelm_exe_name} PRIVATE SOOT_MODEL)
    target_compile_definitions(${pelelm_exe_name} PRIVATE NUM_SOOT_MOMENTS=${PELEMP_NUM_SOOT_MOMENTS})
    target_compile_definitions(${pelelm_exe_name} PRIVATE SOOT_PELE_LM)
    set(SOOT_MOMENTS_VALUES 3 6)
    if(NOT PELEMP_NUM_SOOT_MOMENTS IN_LIST SOOT_MOMENTS_VALUES)
      message(FATAL_ERROR "NUM_SOOT_MOMENTS must be either 3 or 6")
    endif()
    target_sources(${pelelm_exe_name} PRIVATE
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_react.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_derive.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/Constants_Soot.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootData.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootReactions.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_derive.H)
    target_include_directories(${pelelm_exe_name} PRIVATE ${PELEMP_SRC_DIR}/Soot_Models)
  endif()
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
  set_source_files_properties(${PELELM_MECHANISM_DIR}/mechanism.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELELM_MECHANISM_DIR})
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Evaluation)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Utility)
  target_include_directories(${pelelm_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Mechanism/Models)

  target_sources(${pelelm_exe_name}
    PRIVATE
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.cpp
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorTypes.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.H
      ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.cpp
  )

  target_include_directories(${pelelm_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions)
  target_link_libraries(${pelelm_exe_name} PRIVATE sundials_arkode sundials_cvode)

  if(PELELM_ENABLE_CUDA OR PELELM_ENABLE_HIP)
    target_sources(${pelelm_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/AMReX_SUNMemory.cpp
                                              ${PELE_PHYSICS_SRC_DIR}/Reactions/AMReX_SUNMemory.H)
  endif()

  if(PELELM_ENABLE_CUDA)
    target_link_libraries(${pelelm_exe_name} PRIVATE sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
  elseif(PELELM_ENABLE_HIP)
    target_link_libraries(${pelelm_exe_name} PRIVATE sundials_nvechip)
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
       ${SRC_DIR}/Particle.cpp
       ${SRC_DIR}/Soot.cpp
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
       ${IAMR_SRC_DIR}/PROJECTION_F.H
       ${IAMR_SRC_DIR}/NAVIERSTOKES_F.H
       ${IAMR_SRC_DIR}/NS_util.cpp
       ${IAMR_SRC_DIR}/NS_util.H
       ${IAMR_SRC_DIR}/Src_${PELELM_DIM}d/NAVIERSTOKES_${PELELM_DIM}D.F90
  )

  target_sources(${pelelm_exe_name}
     PRIVATE
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_plm.cpp
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov.H
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_edge_state_${PELELM_DIM}D.cpp
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_ppm.H
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov.cpp
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_ppm.cpp
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_K.H
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_extrap_vel_to_faces_${PELELM_DIM}D.cpp
       ${AMREX_HYDRO_SRC_DIR}/Godunov/hydro_godunov_plm.H

       ${AMREX_HYDRO_SRC_DIR}/MOL/hydro_mol.H
       ${AMREX_HYDRO_SRC_DIR}/MOL/hydro_mol.cpp
       ${AMREX_HYDRO_SRC_DIR}/MOL/hydro_mol_extrap_vel_to_faces.cpp
       ${AMREX_HYDRO_SRC_DIR}/MOL/hydro_mol_extrap_vel_to_faces_box.cpp
       ${AMREX_HYDRO_SRC_DIR}/MOL/hydro_mol_edge_state_K.H
       ${AMREX_HYDRO_SRC_DIR}/MOL/hydro_mol_edge_state.cpp

       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_utils.H
       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_utils.cpp
       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_extrap_vel_to_faces.cpp
       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_constants.H
       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_compute_fluxes_from_state.cpp
       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_bcs_K.H
       ${AMREX_HYDRO_SRC_DIR}/Utils/hydro_create_umac_grown.cpp

       ${AMREX_HYDRO_SRC_DIR}/Slopes/hydro_slopes_K.H
       ${AMREX_HYDRO_SRC_DIR}/Slopes/hydro_eb_slopes_${PELELM_DIM}D_K.H
  )

  if(PELELM_ENABLE_AMREX_EB)
    target_sources(${pelelm_exe_name}
       PRIVATE
         ${AMREX_HYDRO_SRC_DIR}/Redistribution/hydro_create_itracker_${PELELM_DIM}d.cpp
         ${AMREX_HYDRO_SRC_DIR}/Redistribution/hydro_redistribution.H
         ${AMREX_HYDRO_SRC_DIR}/Redistribution/hydro_redistribution.cpp
         ${AMREX_HYDRO_SRC_DIR}/Redistribution/hydro_state_redistribute.cpp
         ${AMREX_HYDRO_SRC_DIR}/Redistribution/hydro_state_utils.cpp

         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov.H
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_corner_couple.H
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_edge_state_${PELELM_DIM}D.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_extrap_vel_to_faces.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_extrap_vel_to_faces_${PELELM_DIM}D.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_plm.H
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_plm.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_plm_fpu.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBGodunov/hydro_ebgodunov_transverse_${PELELM_DIM}D_K.H

         ${AMREX_HYDRO_SRC_DIR}/EBMOL/hydro_ebmol.H
         ${AMREX_HYDRO_SRC_DIR}/EBMOL/hydro_ebmol.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBMOL/hydro_ebmol_extrap_vel_to_faces.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBMOL/hydro_ebmol_extrap_vel_to_faces_box.cpp
         ${AMREX_HYDRO_SRC_DIR}/EBMOL/hydro_ebmol_edge_state_K.H
         ${AMREX_HYDRO_SRC_DIR}/EBMOL/hydro_ebmol_edge_state.cpp
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
    target_link_libraries(${pelelm_exe_name} PRIVATE $<$<BOOL:${MPI_C_FOUND}>:MPI::MPI_C>)
    target_link_libraries(${pelelm_exe_name} PRIVATE $<$<BOOL:${MPI_Fortran_FOUND}>:MPI::MPI_Fortran>)
  endif()

  #PeleLM include directories
  target_include_directories(${pelelm_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pelelm_exe_name} PRIVATE ${CMAKE_BINARY_DIR})

  #IAMR include directories
  target_include_directories(${pelelm_exe_name} PRIVATE ${IAMR_SRC_DIR})

  #AMReX-Hydro include directories
  target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/Godunov)
  target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/MOL)
  target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/Utils)
  target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/Slopes)
  if(PELELM_ENABLE_AMREX_EB)
     target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/EBMOL)
     target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/EBGodunov)
     target_include_directories(${pelelm_exe_name} PRIVATE ${AMREX_HYDRO_SRC_DIR}/Redistribution)
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
