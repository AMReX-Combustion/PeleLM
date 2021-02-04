#=============================================================================
# PeleLM Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

# Set location of gold files according to system/compiler/compiler_version
set(FCOMPARE_GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PeleLMGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})

if(PELELM_ENABLE_FCOMPARE_FOR_TESTS)
  message(STATUS "Test golds directory for fcompare: ${FCOMPARE_GOLD_FILES_DIRECTORY}")
endif()

# Have CMake discover the number of cores on the node
#include(ProcessorCount)
#ProcessorCount(PROCESSES)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

macro(setup_test)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/PeleLM-${TEST_EXE_DIR})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_EXE_DIR}/tests/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Find fcompare
    if(PELELM_ENABLE_FCOMPARE_FOR_TESTS)
      set(FCOMPARE ${CMAKE_BINARY_DIR}/Submodules/amrex/Tools/Plotfile/fcompare)
    endif()
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    # Copy files to test working directory
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests
    set(RUNTIME_OPTIONS "amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1")
    if(PELELM_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(PELELM_ENABLE_MPI)
      if(PELELM_ENABLE_CUDA)
        set(PELELM_NP 2) # 1 rank per GPU on Eagle
      else()
        set(PELELM_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELELM_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(PELELM_NP 1)
      unset(MPI_COMMANDS)
    endif()
    # Use fcompare to test diffs in plots against gold files
    if(PELELM_ENABLE_FCOMPARE_FOR_TESTS)
      set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUNTIME_OPTIONS "max_step=10 ${RUNTIME_OPTIONS}")
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${FCOMPARE_COMMAND}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 14400 PROCESSORS ${PELELM_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME TEST_EXE_DIR)
    add_test_r(${TEST_NAME} ${TEST_EXE_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no-ci")
endfunction(add_test_re)

# Standard unit test
function(add_test_u TEST_NAME)
    setup_test()
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/UnitTests/PeleLM-UnitTests)
    set(PELELM_NP 1)
    if(PELELM_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELELM_NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 1800 PROCESSORS ${PELELM_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "unit")
endfunction(add_test_u)

#=============================================================================
# Regression tests
#=============================================================================

# Run in CI
if(PELELM_DIM GREATER 1)
  add_test_r(flamesheet-1 FlameSheet)
  add_test_r(eb-flamepastcylinder-1 EB_FlamePastCylinder)
  #if(PELELM_ENABLE_AMREX_EB)
  #  add_test_r(eb-c4 EB-C4-5)
  #endif()
endif()

# Not run in CI
#if(PELELM_DIM GREATER 1)
#  add_test_re(tg-3 TG)
#  if(PELELM_ENABLE_AMREX_EB)
#    add_test_re(eb-c3 EB-C3)
#  endif()
#endif()

#=============================================================================
# Verification tests
#=============================================================================

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit-tests)

#=============================================================================
# Performance tests
#=============================================================================
