# Local source and header files
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(INC_DIR ${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE SRC_FILES ${SRC_DIR}/*.cpp)
list(REMOVE_ITEM SRC_FILES ${SRC_DIR}/main.cpp)
file(GLOB_RECURSE HEADER_FILES ${INC_DIR}/*.hpp)

# Target definition
blt_add_library(
    NAME quandary_lib
    SOURCES ${SRC_FILES}
    HEADERS ${HEADER_FILES}
    INCLUDES ${INC_DIR} ${CMAKE_BINARY_DIR}/include
    DEPENDS_ON blt::mpi tomlplusplus::tomlplusplus
)

blt_add_executable(
    NAME quandary
    SOURCES ${SRC_DIR}/main.cpp
    DEPENDS_ON quandary_lib
)

# Toml settings
target_compile_definitions(quandary_lib PUBLIC TOML_ENABLE_UNRELEASED_FEATURES=1)

# PETSc
find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
target_link_libraries(quandary_lib PUBLIC PkgConfig::PETSC)

# Put executable in root directory
set_target_properties(
    quandary
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}
)

# Optional: Sanity check
option(SANITY_CHECK "Run sanity checks" OFF)
if(SANITY_CHECK)
    target_compile_options(quandary_lib PUBLIC -DSANITY_CHECK)
endif()

# Optional: Link to SLEPC
option(WITH_SLEPC "Link to SLEPC" OFF)
if(WITH_SLEPC)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(SLEPC REQUIRED IMPORTED_TARGET SLEPc)
    target_link_libraries(quandary_lib PUBLIC PkgConfig::SLEPC)
    target_compile_options(quandary_lib PUBLIC -DWITH_SLEPC)
endif()

install(
    TARGETS quandary
    RUNTIME DESTINATION bin
)
