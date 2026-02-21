set(QUANDARY_PYTHON_INSTALL_DIR "quandary" CACHE PATH "Install destination for Python module")

nanobind_add_module(_quandary_impl STABLE_ABI bindings.cpp)

target_link_libraries(_quandary_impl PRIVATE quandary_lib)

# Generate type stubs for IDE autocomplete
nanobind_add_stub(
    _quandary_impl_stub
    MODULE _quandary_impl
    OUTPUT _quandary_impl.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:_quandary_impl>
    DEPENDS _quandary_impl
    MARKER_FILE py.typed
)

install(TARGETS _quandary_impl LIBRARY DESTINATION ${QUANDARY_PYTHON_INSTALL_DIR})
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/py.typed
    ${CMAKE_CURRENT_BINARY_DIR}/_quandary_impl.pyi
    DESTINATION ${QUANDARY_PYTHON_INSTALL_DIR}
)
