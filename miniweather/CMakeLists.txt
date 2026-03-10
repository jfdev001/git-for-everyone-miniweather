cmake_minimum_required(VERSION 3.0)
project(miniWeather Fortran)

enable_testing()
option(BUILD_TESTING "Build tests" OFF)

include(utils.cmake)


############################################################
## Set Parameters
############################################################
if ("${NX}" STREQUAL "")
  SET(NX 100)
endif()
if ("${NZ}" STREQUAL "")
  SET(NZ 50)
endif()
if ("${SIM_TIME}" STREQUAL "")
  SET(SIM_TIME 1000)
endif()
if ("${OUT_FREQ}" STREQUAL "")
  SET(OUT_FREQ 10)
endif()
if ("${DATA_SPEC}" STREQUAL "")
  SET(DATA_SPEC DATA_SPEC_THERMAL)
endif()
SET(EXE_DEFS "-D_NX=${NX} -D_NZ=${NZ} -D_SIM_TIME=${SIM_TIME} -D_OUT_FREQ=${OUT_FREQ} -D_DATA_SPEC=${DATA_SPEC}")
SET(TEST_DEFS "-D_NX=100 -D_NZ=50 -D_SIM_TIME=400 -D_OUT_FREQ=400 -D_DATA_SPEC=DATA_SPEC_THERMAL")
if ("${FFLAGS}" MATCHES ".*SINGLE_PREC.*")
  set(tol "1.e-5")
  message(STATUS "Using single precision. Test tolerance ${tol}")
else()
  set(tol "1.e-13")
  message(STATUS "Using double precision. Test tolerance ${tol}")
endif()


############################################################
## Append FFLAGS
############################################################
SET(CMAKE_Fortran_FLAGS "${FFLAGS}")


############################################################
## Compile the serial version
############################################################
add_executable(serial miniWeather_serial.F90)
set_target_properties(serial PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")
if(${BUILD_TESTING})
  add_executable(serial_test miniWeather_serial.F90) 
  set_target_properties(serial_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS}")
endif()

if (NOT ("${LDFLAGS}" STREQUAL "") )
  target_link_libraries(serial      "${LDFLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(serial_test "${LDFLAGS}")
  endif()
endif()

if (NOT ("${SERIAL_LINK_FLAGS}" STREQUAL "") )
  target_link_libraries(serial      "${SERIAL_LINK_FLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(serial_test "${SERIAL_LINK_FLAGS}")
  endif()
endif()

if(${BUILD_TESTING})
    add_test(NAME SERIAL_TEST COMMAND ./check_output.sh ./serial_test ${tol} 4.5e-5 ) 
endif()


############################################################
## Compile the MPI version
############################################################
add_executable(mpi miniWeather_mpi.F90)
set_target_properties(mpi PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")

if(${BUILD_TESTING})
    add_executable(mpi_test miniWeather_mpi.F90)
    set_target_properties(mpi_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS}")
endif()

if (NOT ("${LDFLAGS}" STREQUAL "") )
  target_link_libraries(mpi      "${LDFLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(mpi_test "${LDFLAGS}")
  endif()
endif()

if (NOT ("${MPI_LINK_FLAGS}" STREQUAL "") )
  target_link_libraries(mpi      "${MPI_LINK_FLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(mpi_test "${MPI_LINK_FLAGS}")
  endif()
endif()

if(${BUILD_TESTING})
  add_test(NAME MPI_TEST COMMAND ./check_output.sh ./mpi_test ${tol} 4.5e-5 ) 
endif()


############################################################
## Compile the MPI + OpenMP version
############################################################
if (NOT ("${OPENMP_FLAGS}" STREQUAL "") )
  add_executable(openmp miniWeather_mpi_openmp.F90)
  set_target_properties(openmp 
    PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${OPENMP_FLAGS}")

  if(${BUILD_TESTING})
    add_executable(openmp_test miniWeather_mpi_openmp.F90)
    set_target_properties(openmp_test 
      PROPERTIES COMPILE_FLAGS "${TEST_DEFS} ${OPENMP_FLAGS}")
  endif()

  if (NOT ("${LDFLAGS}" STREQUAL "") )
    target_link_libraries(openmp      "${LDFLAGS}")
    if(${BUILD_TESTING})
      target_link_libraries(openmp_test "${LDFLAGS}")
    endif()
  endif()

  if ("${OPENMP_LINK_FLAGS}" STREQUAL "")
    SET(OPENMP_LINK_FLAGS ${OPENMP_FLAGS})
  endif()

  target_link_libraries(openmp      "${OPENMP_LINK_FLAGS}")

  if(${BUILD_TESTING})
    target_link_libraries(openmp_test "${OPENMP_LINK_FLAGS}")
    add_test(NAME OPENMP_TEST COMMAND ./check_output.sh ./openmp_test ${tol} 4.5e-5 ) 
  endif()

endif()


############################################################
## Compile the MPI + OpenACC version
############################################################
if (NOT ("${OPENACC_FLAGS}" STREQUAL "") )
  add_executable(openacc miniWeather_mpi_openacc.F90)
  set_target_properties(openacc PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${OPENACC_FLAGS}")

  add_executable(openacc_test miniWeather_mpi_openacc.F90)
  set_target_properties(openacc_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS} ${OPENACC_FLAGS}")

  if (NOT ("${LDFLAGS}" STREQUAL "") )
    target_link_libraries(openacc      "${LDFLAGS}")
    target_link_libraries(openacc_test "${LDFLAGS}")
  endif()
  if ("${OPENACC_LINK_FLAGS}" STREQUAL "")
    SET(OPENACC_LINK_FLAGS ${OPENACC_FLAGS})
  endif()
  target_link_libraries(openacc      "${OPENACC_LINK_FLAGS}")
  target_link_libraries(openacc_test "${OPENACC_LINK_FLAGS}")

  add_test(NAME OPENACC_TEST COMMAND ./check_output.sh ./openacc_test ${tol} 4.5e-5 ) 
endif()
