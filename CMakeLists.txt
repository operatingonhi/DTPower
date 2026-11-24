cmake_minimum_required(VERSION 3.10)                                              
project(EngineSimulation)
set(CMAKE_CXX_STANDARD 17)                                                                                                                                          
# Detect if we're cross-compiling for Windows 
if(CMAKE_CROSSCOMPILING AND WIN32)
	set(CROSS_COMPILING_WINDOWS TRUE)
	message(STATUS "Cross-compiling for Windows")
else()
	set(CROSS_COMPILING_WINDOWS FALSE)
endif()                                                                                                                                                             
# Add debug and profiling options
option(ENABLE_PROFILING "Enable detailed profiling output" OFF)
option(ENABLE_DEBUG "Enable debug output and checks" OFF)   

if(ENABLE_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)                                                                                                                                                                                                                               
    add_compile_definitions(DEBUG_MODE)                                                                                                                                                                                                                       
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")                                                                                                                                                                                                         
endif()                                                                                                                                                                                                                                                       

if(ENABLE_PROFILING)
    add_compile_definitions(PROFILING_ENABLED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg") 
endif()

                                                                                  
if(CROSS_COMPILING_WINDOWS)
    # Use pre-compiled Cantera for Windows cross-compilation
    #set(CANTERA_ROOT "/home/dev/workspace/cantera")
    #set(CANTERA_INCLUDE_DIRS "${CANTERA_ROOT}/include")
    #set(CANTERA_LIBRARIES "${CANTERA_ROOT}/build/lib/libcantera.a")
    
    set(YAML_CPP_ROOT "/home/dev/workspace/yaml-cpp")
    set(YAML_CPP_INCLUDE_DIRS "${YAML_CPP_ROOT}/include/")
    set(YAML_CPP_LIBRARIES "${YAML_CPP_ROOT}/libyaml-cpp.a")
    
    # Windows-specific libraries that Cantera might need
    set(PLATFORM_LIBRARIES ws2_32 advapi32)
    
    # Static linking flags for MinGW
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static")
    
    # Workaround for missing mutex includes in Cantera headers
    add_compile_definitions(_WIN32_WINNT=0x0601)
    add_compile_definitions(YAML_CPP_STATIC_DEFINE)
    add_compile_options(-include mutex)
    
else()
    # Native compilation - use pkg-config
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CANTERA REQUIRED cantera)

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    
    pkg_check_modules(YAML_CPP REQUIRED yaml-cpp)
    set(PLATFORM_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

include_directories(include)                                                      
                                                                                  
set(SOURCES                                                                       
    src/main.cpp                                                                  
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
                                                                                  
add_executable(engine_sim ${SOURCES})                                             
                                                                                  
target_link_directories(engine_sim PRIVATE
    ${CANTERA_LIBRARY_DIRS}
)

target_link_libraries(engine_sim 
    ${CANTERA_LIBRARIES} 
    ${YAML_CPP_LIBRARIES}
    ${PLATFORM_LIBRARIES}
)
target_include_directories(engine_sim PRIVATE 
    ${CANTERA_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIRS}
)
if(CROSS_COMPILING_WINDOWS)
    target_compile_definitions(engine_sim PRIVATE 
        _WIN32_WINNT=0x0601
        YAML_CPP_STATIC_DEFINE
    )
    target_compile_options(engine_sim PRIVATE -include mutex)
else()
    target_compile_options(engine_sim PRIVATE 
        ${CANTERA_CFLAGS_OTHER}
        ${YAML_CPP_CFLAGS_OTHER}
    )
endif()

# Physical regression checks for the equilibrium thermodynamic updates.
include(CTest)
if(BUILD_TESTING AND NOT CROSS_COMPILING_WINDOWS)
    add_executable(equilibrium_expansion_test
        tests/equilibrium_expansion.cpp
        src/Cylinder.cpp src/CanteraInterface.cpp src/WoschniHeatTransfer.cpp)
    target_include_directories(equilibrium_expansion_test PRIVATE
        ${CANTERA_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIRS})
    target_link_directories(equilibrium_expansion_test PRIVATE ${CANTERA_LIBRARY_DIRS})
    target_link_libraries(equilibrium_expansion_test PRIVATE
        ${CANTERA_LIBRARIES} ${YAML_CPP_LIBRARIES} ${PLATFORM_LIBRARIES})
    target_compile_options(equilibrium_expansion_test PRIVATE
        ${CANTERA_CFLAGS_OTHER} ${YAML_CPP_CFLAGS_OTHER})
    add_test(NAME equilibrium_expansion COMMAND equilibrium_expansion_test)
endif()
