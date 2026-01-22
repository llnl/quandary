nanobind_add_module(quandary_ext STABLE_ABI bindings.cpp)

target_link_libraries(quandary_ext PRIVATE quandary_lib)

install(TARGETS quandary_ext LIBRARY DESTINATION quandary)
