#include <iostream>
#include "cxxopts.hpp"

using namespace std;

void cif2ztb(string in_file, string out_file, 
             float rcut, float k_ewald, float spacing, bool use_fractional_basis, bool out_binary);

cxxopts::ParseResult parse(int argc, char* argv[]) {
  try {
    cxxopts::Options options(argv[0], "Calculates tabulated potentials of material structures");
    options
      .positional_help("[optional args]")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
      ("i,input", "Input file name", cxxopts::value<std::string>())
      ("o,output", "Output file name", cxxopts::value<std::string>()->default_value("zeolite.ztbx"))
      ("s,spacing", "Tabulation spacing in angstroms", cxxopts::value<float>()->default_value("0.1"))
      ("r,rcut", "Cutoff distance in angstroms", cxxopts::value<float>()->default_value("14.0"))
      ("k,kewald", "Ewald parameter, 0 for no Ewald summation, -1 for automatic", cxxopts::value<float>()->default_value("-1"))
      ("f,fractional_basis", "Create the tabulated potential using directions in fractional coordinates",
        cxxopts::value<bool>()->default_value("false"))
      ("b,output_binary", "Outputs a binary file instead of CSV to be read by MCCCS-MN, overrides -f to be True",
        cxxopts::value<bool>()->default_value("false"))
      ("h,help", "Print help")
    ;

    options.parse_positional({"input"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      cout << options.help({"", "Group"}) << endl;
      exit(0);
    }

    return result;

  } catch (const cxxopts::OptionException& e) {
    cout << "error parsing options: " << e.what() << endl;
    exit(1);
  }
}

int main(int argc, char* argv[]) {
  auto result = parse(argc, argv);
  auto arguments = result.arguments();

  string in_file = result["input"].as<string>();
  string out_file = result["output"].as<string>();
  float rcut, kewald, spacing;
  bool frac_basis;
  rcut = result["rcut"].as<float>();
  spacing = result["spacing"].as<float>();
  kewald = result["kewald"].as<float>();
  if (kewald == -1) kewald = 3.2 / rcut;
  
  cif2ztb(in_file, out_file, rcut, kewald, spacing, result["fractional_basis"].as<bool>(), result["output_binary"].as<bool>());
  cout << "Finished potential energy tabulation." << endl;

  return 0;
}