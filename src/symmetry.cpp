#include <sstream>
#include "symmetry.h"

void SymmetryOperation::read_operation_expr(string str, int axis) {
  bool is_denominator = false;
  float current_bias = 0;
  int current_weight = 0;
  for (char const &c: str) {
    switch (c) {
      case '+':
        current_weight = 1;
        break;
      case '-': 
        current_weight = -1;
        break;
      case 'x':
        if (current_weight == 0) current_weight = 1;
        weight[3 * axis] = current_weight;
        break;
      case 'y':
        if (current_weight == 0) current_weight = 1;
        weight[3 * axis + 1] = current_weight;
        break;
      case 'z':
        if (current_weight == 0) current_weight = 1;
        weight[3 * axis + 2] = current_weight;
        break;
      case '/':
        is_denominator = true;
        break;
      case ' ':
        break;
      default:
        if (c >= '1' && c <= '9') {
          if (is_denominator)
            current_bias /= (int) c - '0';
          else 
            current_bias = (int) c - '0';
        } else {
          //throw "Invalid symmetry operation format!";
        }
    }
  }
  bias[axis] = current_bias;
}

void SymmetryOperation::set_operation(string line) {
  bool quoted = false;
  stringstream ss;
  ss << line;
  if (ss.peek() == '\'')
    quoted = true;
  string substr;
  getline(ss, substr, ',' );
  if (quoted)
    substr.erase(0, 1);
  read_operation_expr(substr, 0);
  getline(ss, substr, ',' );
  read_operation_expr(substr, 1);
  getline(ss, substr, ',' );
  if (quoted)
    substr.erase(substr.end() - 1);
  read_operation_expr(substr, 2);
}

vec3 SymmetryOperation::apply(vec3 coords) {
  vec3 new_coords = vec3(bias);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      new_coords[i] += weight[3 * i + j] * coords[j];
    }
    if (new_coords[i] < 0) new_coords[i]++;
    if (new_coords[i] >= 1) new_coords[i]--;
  }
  return new_coords;
} 