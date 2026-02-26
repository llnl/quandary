set(QUANDARY_PYTHON_INSTALL_DIR "quandary" CACHE PATH "Install destination for Python module")

# STABLE_ABI: extension works across Python minor versions without recompilation
nanobind_add_module(_quandary_impl STABLE_ABI bindings.cpp)

target_link_libraries(_quandary_impl PRIVATE quandary_lib)

set_target_properties(_quandary_impl PROPERTIES
    SKIP_BUILD_RPATH FALSE
    BUILD_WITH_INSTALL_RPATH FALSE
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

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
