#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "unitcell.h"

void UnitCell::read_from_cif(string filename) {
  string line, token;
  ifstream file(filename);

  int loop_state = 0;
  int loop_col_index = 0;
  int atom_type_col;
  int axis_cols[3] = {-1, -1, -1};
  if (file.is_open()) {
    while (getline(file, line)) {
      if (line[0] == '#' || line.empty()) continue;
      stringstream ss;
      ss << line;
      ss >> token;
      if (token == "loop_") {
        loop_state = 1;
        continue;
      }
      if (token[0] == '_') {
        if (loop_state == 2)
          loop_state = 0;
        if (token == "_cell_length_a") 
          ss >> cell_length[0];
        else if (token == "_cell_length_b")
          ss >> cell_length[1];
        else if (token == "_cell_length_c")
          ss >> cell_length[2];
        else if (token == "_cell_angle_alpha")
          ss >> cell_angle[0];
        else if (token == "_cell_angle_beta") 
          ss >> cell_angle[1];
        else if (token == "_cell_angle_gamma") 
          ss >> cell_angle[2];
        else if (token == "_space_group_symop_operation_xyz" 
              || token == "_symmetry_equiv_pos_as_xyz") {
          symmetry = get_symmetry(file);
          loop_state = 2;
          loop_col_index = 0;
        }
        
        if (token == "_atom_site_type_symbol")
          atom_type_col = loop_col_index;
          //cout <<atom_type_col<< endl;
        if (token == "_atom_site_fract_x")
          axis_cols[0] = loop_col_index;
        if (token == "_atom_site_fract_y")
          axis_cols[1] = loop_col_index;
        if (token == "_atom_site_fract_z")
          axis_cols[2] = loop_col_index;
        if (loop_state == 1) {
          loop_col_index++;
          //cout <<loop_col_index<<' ' << token << endl;
        }
      } else {
        if (loop_state == 1)
          loop_state = 2;
        if (loop_state == 2) {
          string str = token;
          int i_col = 0;
          bool has_next = true;
          vec3 coords;
          while(has_next) {
            if (i_col == atom_type_col) {
              int id = 0;
              while (id < atomtypes.size() && atomtypes[id] != str)
                id++;
              if (id == atomtypes.size())
                atomtypes.push_back(str);
              reduced_atom_index.push_back(id);
            }
            for (int ax = 0; ax < 3; ax++) {
              if (i_col == axis_cols[ax])
                coords[ax] = stof(str);
            }
            i_col++;
            if (ss >> str) 
              has_next = true;
            else
              has_next = false;
          }
          reduced_coordinates.push_back(coords);
        }
      }   
    }
    //cout << cell_length[0] << '\t' << cell_length[1] << '\t' << cell_length[2] << '\n';
    //cout << cell_angle[0] << '\t' << cell_angle[1] << '\t' << cell_angle[2] << '\n';
    //cout << axis_cols[0] << '\t' << axis_cols[1] << '\t' << axis_cols[2] << '\n';
    file.close();
    calculate_all_coordinates();
  }
  else cout << "Unable to open file"; 
}

void UnitCell::get_reduced_coordinates(vector<vec3> *dest_coordinates, 
          vector<int> *dest_atoms) {
  *dest_coordinates = reduced_coordinates;
  *dest_atoms = reduced_atom_index;
}

void UnitCell::get_all_coordinates(vector<vec3> *dest_coordinates, 
          vector<int> *dest_atoms, bool use_frac) {
  *dest_coordinates = coordinates;
  *dest_atoms = atom_index;
  if (!use_frac)
    for (int i = 0; i < dest_coordinates->size(); i++) {
      (*dest_coordinates)[i] = frac_to_abs((*dest_coordinates)[i] * cell_length);
    }
}

vec3 UnitCell::frac_to_abs(vec3 frac_vector) {
  return frac_to_abs_static(frac_vector, cell_angle);
}

vec3 frac_to_abs_static(vec3 frac_vector, vec3 cell_angle) {
  float cos_alpha, cos_beta, cos_gamma, sin_gamma, c_y, c_z;
  cos_alpha = cos(RAD(cell_angle[0]));
  cos_beta = cos(RAD(cell_angle[1]));
  cos_gamma = cos(RAD(cell_angle[2]));
  sin_gamma = sin(RAD(cell_angle[2]));
  c_y = frac_vector[2] * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
  c_z = frac_vector[2] / sin_gamma * sqrt(
          1 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
          + 2 * cos_alpha * cos_beta * cos_gamma
        );
  vec3 abs_vector;
  abs_vector[0] = frac_vector[0] + frac_vector[1] * cos_gamma 
                                 + frac_vector[2] * cos_beta;
  abs_vector[1] = frac_vector[1] * sin_gamma + c_y;
  abs_vector[2] = c_z;
  return abs_vector;
}

// Transformation matirx in upper triangular form, 6 non-zero elements
float* UnitCell::transformation_matrix() {
  float cos_alpha, cos_beta, cos_gamma, sin_gamma, c_y, c_z;
  float* trig_mat = new float[9];
  for (int i = 0; i <= 9; i++) trig_mat[i] = 0;
  cos_alpha = cos(RAD(cell_angle[0]));
  cos_beta = cos(RAD(cell_angle[1]));
  cos_gamma = cos(RAD(cell_angle[2]));
  sin_gamma = sin(RAD(cell_angle[2]));
  c_y = cell_length[2] * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
  c_z = cell_length[2] / sin_gamma * sqrt(
          1 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
          + 2 * cos_alpha * cos_beta * cos_gamma
        );
  trig_mat[0] = cell_length[0];
  trig_mat[1] = cell_length[1] * cos_gamma;
  trig_mat[2] = cell_length[2] * cos_beta;
  trig_mat[4] = cell_length[1] * sin_gamma;
  trig_mat[5] = c_y;
  trig_mat[8] = c_z;
  return trig_mat;
}


// Invert an upper triangular matrix: [[a,b,c],[0,d,e],[0,0,f]]
float* invert_upper_trig_matrix(float* mat) {
  float* inv_mat = new float[9];
  for (int i = 0; i <= 9; i++) inv_mat[i] = 0;
  inv_mat[0] = 1 / mat[0];
  inv_mat[1] = -mat[1] / mat[0] / mat[4];
  inv_mat[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / (mat[0] * mat[4] * mat[8]);
  inv_mat[4] = 1 / mat[4];
  inv_mat[5] = -mat[5] / mat[4] / mat[8];
  inv_mat[8] = 1 / mat[8];
  return inv_mat;
}


void UnitCell::calculate_all_coordinates() {
  const float TOL = 0.0001;
  coordinates = reduced_coordinates;
  atom_index = reduced_atom_index;
  for (int i = 0; i < reduced_atom_index.size(); i++) {
    vec3 cur_coords = reduced_coordinates[i];
    vector<vec3> cache;
    cache.push_back(cur_coords);
    for (SymmetryOperation symm : symmetry) {
      bool is_new_position = true;
      vec3 new_coords = symm.apply(cur_coords); 
      for (vec3 existing_coords : cache) {
         vec3 diff = new_coords - existing_coords;
         if (diff.squared_length() < TOL * TOL) {
          is_new_position = false;
          break;
        }
      }
      if (is_new_position) {
        coordinates.push_back(new_coords);
        atom_index.push_back(reduced_atom_index[i]);
        cache.push_back(new_coords);
      }
    }
  }
}

// Do minimum image beforehand, do not fold coordinates
// Does not work for ewald sum!
void UnitCell::get_mimage_coordinates(vector<vec3> *dest_coordinates, 
          vector<int> *dest_atoms, float rcut) {
  int nx, ny, nz, size;
  vec3 ex, ey, ez;
  get_all_coordinates(dest_coordinates, dest_atoms, false);  
  nx = rcut / cell_length[0] + 1;
  ny = rcut / cell_length[1] + 1;
  nz = rcut / cell_length[2] + 1;
  ex = frac_to_abs(vec3(cell_length[0], 0, 0));
  ey = frac_to_abs(vec3(0, cell_length[1], 0));
  ez = frac_to_abs(vec3(0, 0, cell_length[2]));
  size = dest_atoms->size();
  for (int x = 0; x < size; x++) 
    for (int i = -nx; i <= nx; i++)
      for (int j = -ny; j <= ny; j++)
        for (int k = -nz; k <= nz; k++)
          if (i != 0 || j != 0 || k != 0) {
            dest_coordinates->push_back((*dest_coordinates)[x] + i * ex + j * ey + k * ez);
            dest_atoms->push_back((*dest_atoms)[x]);
          }
}

float UnitCell::volume() {
  vec3 ex, ey, ez;
  ex = frac_to_abs(vec3(cell_length[0], 0, 0));
  ey = frac_to_abs(vec3(0, cell_length[1], 0));
  ez = frac_to_abs(vec3(0, 0, cell_length[2]));
  return dot(ex, cross(ey, ez));
}

vector<SymmetryOperation> get_symmetry(ifstream &file) {
  string line;
  vector<SymmetryOperation> symm_ops;
  while (file.peek() != '_' && file.peek() != 'l') {
    getline(file, line);
    if (line.empty()) continue;
    //cout << line << endl;
    SymmetryOperation cur_symm_op;
    cur_symm_op.set_operation(line);
    symm_ops.push_back(cur_symm_op);
  }
  return symm_ops;
}

UnitCell UnitCell::duplicate(vec3 dup) {
  int size;
  UnitCell cell_duplicate;
  cell_duplicate.cell_length = dup * cell_length;
  cell_duplicate.cell_angle = cell_angle;
  cell_duplicate.atomtypes = atomtypes;
  size = atom_index.size();
  for (int x = 0; x < size; x++) 
    for (int i = 0; i < dup[0]; i++)
      for (int j = 0; j < dup[1]; j++)
        for (int k = 0; k < dup[2]; k++) {
            cell_duplicate.coordinates.push_back((coordinates[x] + vec3(i, j, k)) / dup);
            cell_duplicate.atom_index.push_back(atom_index[x]);
        }
  return cell_duplicate;
}