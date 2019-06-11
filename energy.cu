#include <iostream>
#include <string>
#include <math.h>
#include "vec3_cuda.h"
#include "unitcell.h"
#define EPS 1e-10

__device__ const float LJSCALE = 2e4;
__device__ const int MAX_ATOMS = 32;



__device__ vec3 get_energy_terms(float r2, float rmin2, float rmax2, float k_ewald) {
    vec3 overlap(1e20, 1e20, 1e20);
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
            tab[2] *= erfcf(k_ewald * r);
    }
    //printf("%f\n", tab[0]);
    return tab;
}

__host__ __device__ 
vec3 frac_to_abs_static_cu(vec3 frac_vector, vec3 cell_angle) {
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

__host__ __device__ 
float* invert_upper_trig_matrix_cu(float* mat) {
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

// returns number of k vectors
__host__ int get_kvectors(vec3* kvectors, float* norms, float k_ewald, float* transform_mat, 
                        float* mat_inv) {
    int n_kvectors = 0;
    int kmin[3] = {0};
    int* kmax = new int[3];
    for (int i = 0; i < 3; i++) {
        kmax[i] = int(transform_mat[4 * i] * k_ewald) + 1;
        if (abs(transform_mat[1]) < EPS || abs(transform_mat[2]) < EPS 
            || abs(transform_mat[5] < EPS))
            kmax[i]++;
        n_kvectors *= 2 * kmax[i] + 1;
    }
    float alpha_2 = k_ewald * k_ewald;
    
    //cudaMallocManaged(&norms, n_kvectors * sizeof(float));
    //cudaMallocManaged(&kvectors, n_kvectors * sizeof(vec3));
    
    float hmax_2 = 4 * PI * PI * alpha_2;
    float * mat_kspace = new float[9];

    for (int i = 0; i < 9; i++) mat_kspace[i] = mat_inv[i] * 2 * PI;
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
                    norms[n_kvectors] = ki_2;
                    kvectors[n_kvectors] = ki;              
                    n_kvectors++;
                    
                }
                
            }
        }
    }
    return n_kvectors;
}

// evald summation from MCCCS-MN; returns electrostatic energy times unit cell volume
__device__ void ewald_sum(float* v_kspace, vec3 r, vec3* atoms_pos, int* atoms_label, int natoms, int atom_types, 
                vec3* k_vectors, float* norms, int n_kvectors, float alpha_2) {
    float arg;
    for (int z = 0; z < atom_types; z++) v_kspace[z] = 0;
    for (int i = 0; i < n_kvectors; i++) {
        float sums[MAX_ATOMS];
        for (int z = 0; z < atom_types; z++) sums[z] = 0;
        for (int z = 0; z < natoms; z++) {
            arg = dot(k_vectors[i], r - atoms_pos[z]);
            sums[atoms_label[z]] += cosf(arg);
        }
        for (int z = 0; z < atom_types; z++) {
            v_kspace[z] += sums[z] * expf(-norms[i] / (4 * alpha_2)) / norms[i];
            
        }
        free(sums);
    }
    for (int z = 0; z < atom_types; z++) v_kspace[z] *= 8 * PI;
}

__host__ __device__ uint32_t index_3d(int i, int j, int k, int dims[]) {
    return i * dims[1] * dims[2] + j * dims[2] + k;
}

__host__ __device__ 
void index_1d(int* ijk, uint32_t index_3d, int* dims) {
    ijk[0] = index_3d / (dims[1] * dims[2]);
    index_3d = index_3d % (dims[1] * dims[2]);
    ijk[1] = index_3d / dims[2];
    ijk[2] = index_3d % dims[2];
}

__global__
void create_grid_pos(vec3* grid, float spacing, int* dims, bool use_fractional, 
                          vec3 cell_angle) {
    int ijk[3];
    uint32_t ind_init = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t stride = blockDim.x * gridDim.x;
    for (uint32_t ind = ind_init; ind < dims[0] * dims[1] * dims[2]; ind+=stride) {
        //printf("%d\n", ind);
        index_1d(ijk, ind, dims);
        for (int z = 0; z < 3; z++)
            grid[ind][z] = spacing * ijk[z];
        if (use_fractional)
            grid[ind] = frac_to_abs_static_cu(grid[ind], cell_angle);
    }
}

// Minimum image convention; box length should be larger than 2rcut
__host__ __device__ vec3 fold_mimage(vec3 v, float* transform_mat, float* mat_inv) {
    vec3 v_frac, v_folded;
    //printf("%f,%f,%f\n", v[0], transform_mat[0], mat_inv[0]);
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


__global__
void calculate_grid(vec3* grid, uint32_t gridsize, 
                    vec3* atoms_pos, int* atom_labels, int natoms, int atom_types,
                    float rmin, float rmax, bool fold, float* transform_mat, float* mat_inv,
                    float k_ewald, int n_kvectors, vec3* kvectors, float* norms, float cell_vol) {
    float rmin2, rmax2, alpha_2;
    float ewald[MAX_ATOMS];
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    alpha_2 = k_ewald * k_ewald;
    
    // Main loop
    uint32_t ind_init = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t stride = blockDim.x * gridDim.x;
    
    for (uint32_t ind = ind_init; ind < gridsize; ind+=stride) {
        
        vec3 pos = grid[ind];
        for (int i = 0; i < atom_types; i++)
            grid[i * gridsize + ind] = vec3(0, 0, 0);  
        for (int j = 0; j < natoms; j++) {
            vec3 dr = pos - atoms_pos[j];
            if (fold) dr = fold_mimage(dr, transform_mat, mat_inv);
            grid[atom_labels[j] * gridsize + ind] += get_energy_terms(dr.squared_length(), rmin2, rmax2, k_ewald);
            
        } 
        if (k_ewald > 0) {
            ewald_sum(ewald, pos, atoms_pos, atom_labels, natoms, atom_types, 
                            kvectors, norms, n_kvectors, alpha_2);
            for (int i = 0; i < atom_types; i++)
                grid[i * gridsize + ind][2] += ewald[i] / cell_vol;
        }
    }
    free(ewald);
}

