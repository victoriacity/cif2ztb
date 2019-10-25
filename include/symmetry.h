#ifndef SYMMETRYH
#define SYMMETRYH

#include <string>
#include "vec3.h"

using namespace std;

const float CELL_TOL = 1e-5;

class SymmetryOperation {
    int weight[9] = {0};
    float bias[3] = {0.0};
    void read_operation_expr(string str, int axis);

    public:
    void set_operation(string line);
    vec3 apply(vec3 coords);
};

#endif