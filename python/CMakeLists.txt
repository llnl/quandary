set(QUANDARY_PYTHON_INSTALL_DIR "quandary" CACHE PATH "Install destination for Python module")

# STABLE_ABI: extension works across Python minor versions without recompilation
nanobind_add_module(_quandary_impl STABLE_ABI bindings.cpp)

target_link_libraries(_quandary_impl PRIVATE quandary_lib)

set_target_properties(_quandary_impl PROPERTIES
    SKIP_BUILD_RPATH FALSE
    BUILD_WITH_INSTALL_RPATH FALSE
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# PETSc, SLEPc, and MPI come transitively through quandary_lib, but CMake's
# INSTALL_RPATH_USE_LINK_PATH does not capture library directories from
# PkgConfig imported targets or BLT targets that are inherited transitively.
# Explicitly add their library dirs to both BUILD_RPATH and INSTALL_RPATH
# so the nanobind extension can find them at runtime.
if(PETSC_LIBRARY_DIRS)
    set_property(TARGET _quandary_impl APPEND PROPERTY
        BUILD_RPATH "${PETSC_LIBRARY_DIRS}")
    set_property(TARGET _quandary_impl APPEND PROPERTY
        INSTALL_RPATH "${PETSC_LIBRARY_DIRS}")
endif()
if(WITH_SLEPC AND SLEPC_LIBRARY_DIRS)
    set_property(TARGET _quandary_impl APPEND PROPERTY
        BUILD_RPATH "${SLEPC_LIBRARY_DIRS}")
    set_property(TARGET _quandary_impl APPEND PROPERTY
        INSTALL_RPATH "${SLEPC_LIBRARY_DIRS}")
endif()
if(MPI_CXX_LIBRARIES)
    foreach(lib ${MPI_CXX_LIBRARIES})
        if(EXISTS "${lib}")
            get_filename_component(mpi_lib_dir "${lib}" DIRECTORY)
            set_property(TARGET _quandary_impl APPEND PROPERTY
                BUILD_RPATH "${mpi_lib_dir}")
            set_property(TARGET _quandary_impl APPEND PROPERTY
                INSTALL_RPATH "${mpi_lib_dir}")
        endif()
    endforeach()
endif()

# Generate type stubs for IDE autocomplete (requires typing_extensions on Python < 3.11)
option(BUILD_PYTHON_STUBS "Generate Python type stubs (.pyi) for the nanobind module" ON)
if(BUILD_PYTHON_STUBS)
    nanobind_add_stub(
        _quandary_impl_stub
        MODULE _quandary_impl
        OUTPUT _quandary_impl.pyi
        PYTHON_PATH $<TARGET_FILE_DIR:_quandary_impl>
        DEPENDS _quandary_impl
        MARKER_FILE py.typed
    )
endif()

install(TARGETS _quandary_impl LIBRARY DESTINATION ${QUANDARY_PYTHON_INSTALL_DIR})
if(BUILD_PYTHON_STUBS)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/py.typed
        ${CMAKE_CURRENT_BINARY_DIR}/_quandary_impl.pyi
        DESTINATION ${QUANDARY_PYTHON_INSTALL_DIR}
    )
endif()

# Install pure Python source files (only for direct CMake builds, not scikit-build-core)
if(NOT SKBUILD)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/python/quandary/
        DESTINATION ${QUANDARY_PYTHON_INSTALL_DIR}
        FILES_MATCHING PATTERN "*.py"
    )
endif()
