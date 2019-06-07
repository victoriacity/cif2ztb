#ifndef UNITCELLH
#define UNITCELLH

#include <string>
#include <vector>

#include "vec3.h"
#include "symmetry.h"

#define PI 3.14159265
#define RAD(x) (x*PI/180)

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
  void get_all_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms, bool use_frac);
  void get_reduced_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms);
  void get_mimage_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms, float rcut);
  vec3 frac_to_abs(vec3 frac_vector);
  float* transformation_matrix();
  float volume();
};

vector<SymmetryOperation> get_symmetry(ifstream &file);
vec3 frac_to_abs_static(vec3 frac_vector, vec3 cell_angle);
float* invert_upper_trig_matrix(float* mat);
UnitCell duplicate(Unitcell cell, vec3 dup);

#endif

