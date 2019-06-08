#include <iostream>
#include <string>
#include <math.h>
#include "vec3.h"
#include "unitcell.h"

#define EPS 1e-10

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
            tab[2] *= 1 - erf(k_ewald * r);
    }
    return tab;
}

// returns number of k vectors
int get_kvectors(int* kvectors, float* kvector_norms, float* mat_inv, float alpha_2) {
    return 0;
}

// evald summation from MCCCS-MN; returns electrostatic energy times unit cell volume
float* ewald_sum(vec3 r, vec3* atoms_pos, int* atoms_label, int natoms, int atom_types, 
                float* mat_kspace, 
                int kmax[], float alpha_2) {
    int kmin[3] = {0};
    float hmax_2, ki_2, arg;
    float* v_kspace = new float[atom_types];
    for (int z = 0; z < atom_types; z++) v_kspace[z] = 0;
    vec3* tab = new vec3[atom_types];
    vec3 ki;
    hmax_2 = 4 * PI * PI * alpha_2;
    for (int l = 0; l <= kmax[0]; l++) {
        if (l == 0) kmin[1] = 0; else kmin[1] = -kmax[1];
        for (int m = kmin[1]; m <= kmax[1]; m++) {
            if (l == 0 && m == 0) kmin[2] = 1; else kmin[2] = -kmax[2];
            for (int n = kmin[2]; n <= kmax[2]; n++) {
                vec3 ki;
                float ki_2;
                for (int z = 0; z < 3; z++) {
                    ki[z] = l * mat_kspace[z] + m * mat_kspace[3 + z] + n * mat_kspace[6 + z];
                }
                ki_2 = ki.squared_length();
                if (hmax_2 - ki_2 > EPS) {
                    float sums[atom_types] = {0};
                    for (int z = 0; z < natoms; z++) {
                        arg = dot(ki, r - atoms_pos[z]);
                        sums[atoms_label[z]] += cos(arg);
                    }
                    for (int z = 0; z < atom_types; z++) {
                        v_kspace[z] += sums[z] * exp(-ki_2 / (4 * alpha_2)) / ki_2;
                        
                    }
                }
            }
        }
    }
    for (int z = 0; z < atom_types; z++) v_kspace[z] *= 8 * PI;
    return v_kspace;
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
                          vec3 cell_angle) {
    int* ijk = new int[3];
    uint32_t ind;
    ijk = index_1d(0, dims);
    for (ind = 0; ind < dims[0] * dims[1] * dims[2]; ind++) {
        for (int z = 0; z < 3; z++)
            grid[ind][z] = spacing * ijk[z];
        if (use_fractional)
            grid[ind] = frac_to_abs_static(grid[ind], cell_angle);
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

// Minimum image convention; box length should be larger than 2rcut
vec3 fold_mimage(vec3 v, float* transform_mat, float* mat_inv) {
    vec3 v_frac, v_folded;
    for (int i = 0; i < 3; i++) {
        v_frac[i] = mat_inv[3 * i] * v[0] + mat_inv[3 * i + 1] * v[1] + mat_inv[3 * i + 2] * v[2];
        if (v_frac[i] <= -0.5) v_frac[i]++;
        if (v_frac[i] > 0.5) v_frac[i]--;
    }
    for (int i = 0; i < 3; i++)
        v_folded[i] = transform_mat[3 * i] * v_frac[0] 
                        + transform_mat[3 * i + 1] * v_frac[1]
                        + transform_mat[3 * i + 2] * v_frac[2];
    return v_folded;
}


void calculate_grid(vec3** grid, uint32_t gridsize, 
                    vec3* atoms_pos, int* atom_labels, int natoms, int atom_types,
                    float rmin, float rmax, float k_ewald, bool fold,
                    float* transform_mat, float cell_vol) {
    uint32_t ind;
    float rmin2, rmax2, alpha_2;
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    float *mat_inv, *mat_kspace;
    int *kmax;
    mat_inv = invert_upper_trig_matrix(transform_mat);
    if (k_ewald > 0) {
        kmax = new int[3];
        mat_kspace = new float[9];
        for (int i = 0; i < 9; i++) mat_kspace[i] = mat_inv[i] * 2 * PI;
        for (int i = 0; i < 3; i++) {
            kmax[i] = int(transform_mat[4 * i] * k_ewald) + 1;
            if (abs(transform_mat[1]) < EPS || abs(transform_mat[2]) < EPS 
                || abs(transform_mat[5] < EPS))
                kmax[i]++;
        }
        alpha_2 = k_ewald * k_ewald;
    }

    // Main loop
    for (ind = 0; ind < gridsize; ind++) {
        vec3 pos = grid[0][ind];
        for (int i = 0; i < atom_types; i++)
            grid[i][ind] = vec3(0, 0, 0);
        for (int j = 0; j < natoms; j++) {
            vec3 dr = pos - atoms_pos[j];
            if (fold) dr = fold_mimage(dr, transform_mat, mat_inv);
            grid[atom_labels[j]][ind] += get_energy_terms(dr.squared_length(), rmin2, rmax2, k_ewald);      
        }
        if (k_ewald > 0) {
            float* ewald = ewald_sum(pos, atoms_pos, atom_labels, natoms, atom_types, 
                            mat_kspace, kmax, alpha_2);
            for (int i = 0; i < atom_types; i++)
                grid[i][ind][2] += ewald[i] / cell_vol;
        }
        //std::cout << grid[0][ind][2] << endl;
    }
}