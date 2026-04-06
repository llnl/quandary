if(ENABLE_TESTS)
  enable_testing()

  blt_add_executable(
    NAME test_cfgparser
    SOURCES test_cfgparser.cpp
    DEPENDS_ON quandary_lib gtest
    FOLDER tests
  )

  blt_add_test(
    NAME test_cfgparser
    COMMAND test_cfgparser
  )

  blt_add_executable(
    NAME test_tomlparser
    SOURCES test_tomlparser.cpp
    DEPENDS_ON quandary_lib gtest
    FOLDER tests
  )

  blt_add_test(
    NAME test_tomlparser
    COMMAND test_tomlparser
  )
endif()
