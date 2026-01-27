nanobind_add_module(_quandary_impl STABLE_ABI bindings.cpp)

target_link_libraries(_quandary_impl PRIVATE quandary_lib)

# Generate type stubs for IDE autocomplete
nanobind_add_stub(
    _quandary_impl_stub
    MODULE _quandary_impl
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/quandary/_quandary_impl.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:_quandary_impl>
    DEPENDS _quandary_impl
)

install(TARGETS _quandary_impl LIBRARY DESTINATION quandary)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/quandary/_quandary_impl.pyi DESTINATION quandary)
