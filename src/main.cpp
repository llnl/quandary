#include "quandary_core.hpp"
#include "util.hpp"

int main(int argc, char** argv) {
  ParsedArgs args = parseArguments(argc, argv);
  try {
    Config config = Config::fromFile(args.config_filename, args.quietmode);
    return runQuandary(config, args.quietmode, argc, argv, args.petsc_argc, args.petsc_argv.data());
  } catch (const std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
}
