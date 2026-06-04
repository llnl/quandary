#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

#include "util.hpp"

// parseArguments() validates that the config file (argv[1]) exists before
// parsing flags, so tests provide an empty, readable file as the config.
class ArgParserTest : public ::testing::Test {
 protected:
  const std::string config = "argparser_test.cfg";

  void SetUp() override { std::ofstream(config).close(); }
  void TearDown() override { std::remove(config.c_str()); }

  // Calls parseArguments() with argv built from the given tokens.
  ParsedArgs parse(std::vector<std::string> tokens) {
    std::vector<char*> argv;
    for (auto& tok : tokens) argv.push_back(const_cast<char*>(tok.c_str()));
    return parseArguments(static_cast<int>(argv.size()), argv.data());
  }
};

TEST_F(ArgParserTest, ParsesQuietAndQuotedPetscOptions) {
  // The shell strips the quotes the user types and passes the quoted options as
  // a SINGLE argv element containing a space ("-log_view -tao_view"), which
  // parseArguments then splits on whitespace. We reproduce that here directly.
  ParsedArgs args = parse({"quandary", config, "--quiet", "--petsc-options",
                           "-log_view -tao_view"});
  EXPECT_TRUE(args.quietmode);
  ASSERT_EQ(args.petsc_tokens.size(), 2);
  EXPECT_EQ(args.petsc_tokens[0], "-log_view");
  EXPECT_EQ(args.petsc_tokens[1], "-tao_view");
}

// Forgetting the quotes leaves the value-tokens as unrecognized arguments
// (the original crash: PETSc received a value-less -vec_type). This must
// error out rather than silently dropping them.
TEST_F(ArgParserTest, UnquotedPetscOptionsExitWithError) {
  EXPECT_EXIT(parse({"quandary", config, "--petsc-options", "-vec_type",
                     "kokkos", "-mat_type", "aijkokkos"}),
              ::testing::ExitedWithCode(1), "");
}

TEST_F(ArgParserTest, UnknownFlagExitsWithError) {
  EXPECT_EXIT(parse({"quandary", config, "--bogus"}),
              ::testing::ExitedWithCode(1), "");
}

TEST_F(ArgParserTest, PetscOptionsMissingValueExitsWithError) {
  EXPECT_EXIT(parse({"quandary", config, "--petsc-options"}),
              ::testing::ExitedWithCode(1), "");
}

// --help and --version must work in any position (exit 0), not be treated as
// unrecognized arguments when they follow the config file.
TEST_F(ArgParserTest, HelpAfterConfigExitsZero) {
  EXPECT_EXIT(parse({"quandary", config, "--help"}),
              ::testing::ExitedWithCode(0), "");
}

TEST_F(ArgParserTest, VersionAfterConfigExitsZero) {
  EXPECT_EXIT(parse({"quandary", config, "--version"}),
              ::testing::ExitedWithCode(0), "");
}
