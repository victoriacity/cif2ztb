#ifndef UNITCELLH
#define UNITCELLH

#include <string>
#include <vector>

#include "vec3.h"
#include "symmetry.h"

using namespace std;

class UnitCell {
  vector<SymmetryOperation> symmetry;
  vector<vec3> coordinates;
  vector<int> atom_index;

  public:
  vec3 cell_length;
  vec3 cell_angle;
  vector<string> atomtypes;
  void read_from_cif(string filename);
  void get_all_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms);
  void get_reduced_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms);
  vec3 frac_to_abs(vec3 frac_vector);
};

vector<SymmetryOperation> get_symmetry(ifstream &file);

#endif

