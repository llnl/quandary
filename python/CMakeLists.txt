nanobind_add_module(quandary_ext STABLE_ABI bindings.cpp)

target_link_libraries(quandary_ext PRIVATE quandary_lib)

# Generate type stubs for IDE autocomplete
nanobind_add_stub(
    quandary_ext_stub
    MODULE quandary_ext
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/quandary/quandary_ext.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:quandary_ext>
    DEPENDS quandary_ext
)

install(TARGETS quandary_ext LIBRARY DESTINATION quandary)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/quandary/quandary_ext.pyi DESTINATION quandary)
