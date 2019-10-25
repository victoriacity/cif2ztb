#ifndef UNITCELLH
#define UNITCELLH

#include <string>
#include <vector>

#include "vec3.h"
#include "symmetry.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749
#define RAD(x) (x*PI/180)

using namespace std;

class UnitCell {
  vector<SymmetryOperation> symmetry;
  vector<vec3> reduced_coordinates;
  vector<int> reduced_atom_index;
  void calculate_all_coordinates();

  public:
  vec3 cell_length;
  vec3 cell_angle;
  vector<string> atomtypes;

  vector<int> atom_index;
  vector<vec3> coordinates;

  void read_from_cif(string filename);
  void get_all_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms, bool use_frac);
  void get_reduced_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms);
  void get_mimage_coordinates(vector<vec3> *dest_coordinates, vector<int> *dest_atoms, float rcut);
  vec3 frac_to_abs(vec3 frac_vector);
  float* transformation_matrix();
  float volume();
  UnitCell duplicate(vec3 dup);
};

vector<SymmetryOperation> get_symmetry(ifstream &file);
vec3 frac_to_abs_static(vec3 frac_vector, vec3 cell_angle);
float* invert_upper_trig_matrix(float* mat);

#endif

