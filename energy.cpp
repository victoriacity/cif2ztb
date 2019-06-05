#include <iostream>
#include <string>
#include <math.h>
#include "vec3.h"
#include "unitcell.h"

const float LJSCALE = 2e4;
vec3 overlap(1e20, 1e20, 1e20);

vec3 get_energy_terms(float r2, float rmin2, float rmax2, float k_ewald) {
    if (r2 < rmin2) 
        return overlap;
    vec3 tab(0, 0, 0); 
    if (r2 < rmax2) {
        float r = sqrt(r2);
        float r6_inv = LJSCALE / pow(r2, 3);
        tab[0] = r6_inv * r6_inv;
        tab[1] = r6_inv;
        tab[2] = 1 / r;
        if (k_ewald > 0) 
            tab[2] *= erf(k_ewald * r);
    }
    return tab;
}

uint32_t index_3d(int i, int j, int k, int dims[]) {
    return i * dims[1] * dims[2] + j * dims[2] + k;
}

int* index_1d(uint32_t index_3d, int dims[]) {
    int* ijk = new int[3];
    ijk[0] = index_3d / (dims[1] * dims[2]);
    index_3d = index_3d % (dims[1] * dims[2]);
    ijk[1] = index_3d / dims[2];
    ijk[2] = index_3d % dims[2];
    return ijk;
}

void create_grid_pos(vec3* grid, float spacing, int dims[], bool use_fractional, 
                          vec3 (*frac_to_abs)(vec3 frac_vector)) {
    int* ijk = new int[3];
    uint32_t ind;
    ijk = index_1d(0, dims);
    for (ind = 0; ind < dims[0] * dims[1] * dims[2]; ind++) {
        for (int z = 0; z < 2; z++)
            grid[ind][z] = spacing * ijk[z];
        if (use_fractional)
            grid[ind] = frac_to_abs(grid[ind]);
        if ((ind + 1) % (dims[1] * dims[2]) == 0) {
            ijk[0]++;
            ijk[1] = 0;
            ijk[2] = 0;
        } else if ((ind + 1) % dims[2] == 0) {
            ijk[2] = 0;
            ijk[1]++;
        } else        
            ijk[2]++;
    }   
}

void calculate_grid(vec3* grid, uint32_t gridsize, vec3* atoms_pos, int natoms,
                    float rmin, float rmax, float k_ewald) {
    uint32_t ind;
    float rmin2, rmax2;
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    
    for (ind = 0; ind < gridsize; ind++) {
        vec3 pos = grid[ind];
        vec3 tab = vec3(0, 0, 0);
        for (int j = 0; j < natoms; j++) {
            vec3 dr = pos - atoms_pos[j];
            // TODO: minimum image convention
            tab += get_energy_terms(dr.squared_length(), rmin2, rmax2, k_ewald);
        }
        grid[ind] = tab;
    }
}