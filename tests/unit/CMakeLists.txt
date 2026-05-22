if(ENABLE_TESTS)
  enable_testing()

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
