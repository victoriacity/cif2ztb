#include <iostream>
#include <string>
#include "vec3.h"
#include "unitcell.h"

void create_grid_pos(vec3* grid, float spacing, int dims[], bool use_fractional, 
                          vec3 cell_angle);

void calculate_grid(vec3** grid, uint32_t gridsize, 
                    vec3* atoms_pos, int* atom_labels, int natoms, int atom_types,
                    float rmin, float rmax, float k_ewald, bool fold,
                    float* transform_mat, float cell_vol);

int main () {
  vector<vec3> coordinates;
  vector<int> atoms;
  UnitCell cell;
  float rcut = 14.0;
  float k_ewald = 3.2 / rcut;
  float spacing = 0.1;
  int dims[] = {10, 10, 1};
  int atom_types;
  uint32_t natom, size;
  cell.read_from_cif("MFI-0.cif");
  
  cell = cell.duplicate(vec3(2, 2, 3));
  
  cell.get_all_coordinates(&coordinates, &atoms, false);
  //cell.get_mimage_coordinates(&coordinates, &atoms, rcut);
  bool fold = true;
  //for (int i = 0; i < coordinates.size(); i++)
  //  std::cout << atoms[i] << ' ' << coordinates[i] << endl;
  
  atom_types = cell.atomtypes.size();
  natom = coordinates.size();
  size = dims[0] * dims[1] * dims[2];
  vec3** grid = new vec3*[atom_types];
  for (int i = 0; i < atom_types; i++) {
    grid[i] = new vec3[size];
  }
  //k_ewald = 0;
  create_grid_pos(grid[0], spacing, dims, false, cell.cell_angle);
  //for (int i = 0; i < 10; i++) std::cout << grid[0][i] << endl;
  calculate_grid(grid, size, coordinates.data(), atoms.data(), natom, atom_types, 0.05, rcut,
                 k_ewald, fold, cell.transformation_matrix(), cell.volume());  
  for (int i = 0; i < atom_types; i++) {
    for (int j = 0; j < size; j++)
      std::cout << grid[i][j] << ' ';
    std::cout << endl;
  }
  return 0;
}