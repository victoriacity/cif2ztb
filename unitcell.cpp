#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "unitcell.h"

#define PI 3.14159265
#define RAD(x) (x*PI/180)

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
              atom_index.push_back(id);
            }
            for (int ax = 0; ax < 2; ax++) {
              if (i_col == axis_cols[ax])
                coords[ax] = stof(str);
            }
            i_col++;
            if (ss >> str) 
              has_next = true;
            else
              has_next = false;
          }
          coordinates.push_back(coords);
        }
      }   
    }
    //cout << cell_length[0] << '\t' << cell_length[1] << '\t' << cell_length[2] << '\n';
    //cout << cell_angle[0] << '\t' << cell_angle[1] << '\t' << cell_angle[2] << '\n';
    //cout << axis_cols[0] << '\t' << axis_cols[1] << '\t' << axis_cols[2] << '\n';
    file.close();
  }
  else cout << "Unable to open file"; 
}

void UnitCell::get_reduced_coordinates(vector<vec3> *dest_coordinates, 
          vector<int> *dest_atoms) {
  *dest_coordinates = coordinates;
  *dest_atoms = atom_index;
}

vec3 UnitCell::frac_to_abs(vec3 frac_vector) {
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


void UnitCell::get_all_coordinates(vector<vec3> *dest_coordinates, 
          vector<int> *dest_atoms) {
  const float TOL = 0.0001;
  *dest_coordinates = coordinates;
  *dest_atoms = atom_index;
  for (int i = 0; i < atom_index.size(); i++) {
    vec3 cur_coords = coordinates[i];
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
        dest_coordinates->push_back(new_coords);
        dest_atoms->push_back(atom_index[i]);
        cache.push_back(new_coords);
      }
    }
    
  } 
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