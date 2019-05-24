#include <iostream>
#include <string>
#include "vec3.h"
#include "unitcell.h"

int main () {
  vector<vec3> coordinates;
  vector<int> atoms;
  UnitCell cell;
  cell.read_from_cif("MFI-2.cif");
  cell.get_all_coordinates(&coordinates, &atoms);
  cout << atoms.size() << endl;
  cell.get_reduced_coordinates(&coordinates, &atoms);
  cout << atoms.size() << endl;
  for (auto i: cell.atomtypes)
    std::cout << i << ' ';
  cout << endl;
  cout << cell.cell_angle << endl;
  vec3 rel(1, 1, 1);
  vec3 abs = cell.frac_to_abs(rel);
  cout << abs << endl;
  return 0;
}