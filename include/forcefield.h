#ifndef FORCEFIELDH
#define FORCEFIELDH

#include <math.h>
#include <stdlib.h>
#include "vec3.h"

class ForceField {

float repulsion(vec3 r, int atomid);
float attraction(vec3 r, int atomid);
float coulomb(vec3 r, int atomid);
}

#endif