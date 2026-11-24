cmake_minimum_required(VERSION 3.18)
project(EngineSimulation LANGUAGES CXX)

include(GNUInstallDirs)

option(ENABLE_PROFILING "Enable detailed profiling output and gprof instrumentation" OFF)
option(ENABLE_DEBUG "Enable debug output and checks" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(engine_build_options INTERFACE)
target_include_directories(engine_build_options INTERFACE include)

if(ENABLE_DEBUG)
    target_compile_definitions(engine_build_options INTERFACE DEBUG_MODE)
    if(MSVC)
        target_compile_options(engine_build_options INTERFACE /Zi /Od /W4)
        target_link_options(engine_build_options INTERFACE /DEBUG)
    else()
        target_compile_options(engine_build_options INTERFACE -g -O0 -Wall -Wextra)
    endif()
endif()

if(ENABLE_PROFILING)
    target_compile_definitions(engine_build_options INTERFACE PROFILING_ENABLED)
    target_compile_options(engine_build_options INTERFACE -pg)
    target_link_options(engine_build_options INTERFACE -pg)
endif()

if(CMAKE_CROSSCOMPILING AND WIN32)
    # Supply target libraries and headers via the toolchain or CMake cache.
    foreach(variable IN ITEMS CANTERA_INCLUDE_DIRS CANTERA_LIBRARIES
                              YAML_CPP_INCLUDE_DIRS YAML_CPP_LIBRARIES)
        if(NOT ${variable})
            message(FATAL_ERROR "Set ${variable} for the Windows target.")
        endif()
    endforeach()
    target_include_directories(engine_build_options SYSTEM INTERFACE
        ${CANTERA_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIRS})
    if(CANTERA_LIBRARY_DIRS)
        target_link_directories(engine_build_options INTERFACE ${CANTERA_LIBRARY_DIRS})
    endif()
    target_link_libraries(engine_build_options INTERFACE
        ${CANTERA_LIBRARIES} ${YAML_CPP_LIBRARIES} ws2_32 advapi32)
    target_compile_definitions(engine_build_options INTERFACE
        _WIN32_WINNT=0x0601 YAML_CPP_STATIC_DEFINE)
    target_compile_options(engine_build_options INTERFACE -include mutex)
    target_link_options(engine_build_options INTERFACE -static-libgcc -static-libstdc++ -static)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CANTERA REQUIRED IMPORTED_TARGET cantera)
    pkg_check_modules(YAML_CPP REQUIRED IMPORTED_TARGET yaml-cpp)
    find_package(LAPACK REQUIRED)
    target_link_libraries(engine_build_options INTERFACE
        PkgConfig::CANTERA PkgConfig::YAML_CPP LAPACK::LAPACK)
endif()

add_library(engine_core STATIC
    src/Engine.cpp
    src/Cylinder.cpp
    src/WoschniHeatTransfer.cpp
    src/ValveFlow.cpp
    src/FrictionModel.cpp
    src/Crankshaft.cpp
    src/SparkModel.cpp
    src/CanteraInterface.cpp
    src/ConfigLoader.cpp
)
target_link_libraries(engine_core PUBLIC engine_build_options)

add_executable(engine_sim src/main.cpp)
target_link_libraries(engine_sim PRIVATE engine_core)

install(TARGETS engine_sim RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES README.md LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DATADIR}/dtpower/examples
    FILES_MATCHING PATTERN "*_parameters.yaml")

include(CTest)
if(BUILD_TESTING AND NOT (CMAKE_CROSSCOMPILING AND WIN32))
    add_executable(equilibrium_expansion_test tests/equilibrium_expansion.cpp)
    target_link_libraries(equilibrium_expansion_test PRIVATE engine_core)
    add_test(NAME equilibrium_expansion COMMAND equilibrium_expansion_test)

    add_executable(mechanical_models_test tests/mechanical_models.cpp)
    target_link_libraries(mechanical_models_test PRIVATE engine_core)
    add_test(NAME mechanical_models COMMAND mechanical_models_test)

    add_executable(cycle_reuse_test tests/cycle_reuse.cpp)
    target_link_libraries(cycle_reuse_test PRIVATE engine_core)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cycle_reuse)
    add_test(NAME cycle_reuse COMMAND cycle_reuse_test
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/methane-air)
    set_tests_properties(cycle_reuse PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cycle_reuse TIMEOUT 180)

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    foreach(case IN ITEMS convergence operating_range output_interval closed_cylinder invalid_inputs)
        add_test(NAME engine_${case}
            COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/tests/engine_cli.py
                    $<TARGET_FILE:engine_sim> ${CMAKE_CURRENT_SOURCE_DIR}/examples/methane-air ${case})
        set_tests_properties(engine_${case} PROPERTIES TIMEOUT 360)
    endforeach()
endif()
