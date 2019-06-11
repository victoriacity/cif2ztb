#include <iostream>
#include <string>
#include "vec3_cuda.h"
#include "unitcell.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__host__ int get_kvectors(vec3* kvectors, float* norms, float k_ewald, float* transform_mat, 
  float* mat_inv);

__global__
void create_grid_pos(vec3* grid, float spacing, int dims[], bool use_fractional, 
                          vec3 cell_angle);

__global__
void calculate_grid(vec3* grid, uint32_t gridsize, 
          vec3* atoms_pos, int* atom_labels, int natoms, int atom_types,
          float rmin, float rmax, bool fold, float* transform_mat, float* mat_inv,
          float k_ewald, int n_kvectors, vec3* kvectors, float* norms, float cell_vol);

int main () {
  vector<vec3> coordinates;
  vector<int> atoms;
  UnitCell cell;
  float rcut = 14.0;
  float k_ewald = 3.2 / rcut;
  float spacing = 0.1;
  int* dims;
  const uint32_t MAX_KVECTORS = 10000;
  gpuErrchk(cudaMallocManaged(&dims, 3 * sizeof(int)));
  dims[0] = 200;
  dims[1] = 197;
  dims[2] = 131;
  int atom_types;
  uint32_t natom, size;
  cell.read_from_cif("MFI-0.cif");
  //cell.read_from_cif("STF-1.cif");

  cell = cell.duplicate(vec3(2, 2, 3));
  
  cell.get_all_coordinates(&coordinates, &atoms, false);
  bool fold = true;
  //for (int i = 0; i < coordinates.size(); i++)
  //  std::cout << atoms[i] << ' ' << coordinates[i] << endl;
  
  atom_types = cell.atomtypes.size();
  natom = coordinates.size();
  size = dims[0] * dims[1] * dims[2];
  vec3* grid;
  gpuErrchk(cudaMallocManaged(&grid, atom_types * size * sizeof(vec3)));
   
  //for (int i = 0; i < 10; i++) std::cout << grid[i] << endl;
  float* d_mat, *invmat, *d_invmat, *norms;
  vec3* d_coordinates, *kvectors;
  int* d_atoms;
  int n_kvectors;

  gpuErrchk(cudaMalloc(&d_mat, 9 * sizeof(float)));
  gpuErrchk(cudaMalloc(&d_invmat, 9 * sizeof(float)));
  gpuErrchk(cudaMalloc(&d_coordinates, natom * sizeof(vec3)));
  gpuErrchk(cudaMalloc(&d_atoms, natom * sizeof(int)));

  invmat = invert_upper_trig_matrix(cell.transformation_matrix());

  gpuErrchk(cudaMallocManaged(&norms, MAX_KVECTORS * sizeof(float)));
  gpuErrchk(cudaMallocManaged(&kvectors, MAX_KVECTORS * sizeof(vec3)));
  n_kvectors = get_kvectors(kvectors, norms, k_ewald, cell.transformation_matrix(), invmat);

  gpuErrchk(cudaMemcpy(d_mat, cell.transformation_matrix(), 9 * sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_invmat, invmat, 9 * sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_coordinates, coordinates.data(), natom * sizeof(vec3), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_atoms, atoms.data(), natom * sizeof(int), cudaMemcpyHostToDevice));

  uint32_t blockSize = 256;
  uint32_t numBlocks = (size + blockSize - 1) / blockSize;
  //k_ewald = 0;
  create_grid_pos<<<numBlocks, blockSize>>>(grid, spacing, dims, true, cell.cell_angle);

  calculate_grid<<<numBlocks, blockSize>>>(grid, size, d_coordinates, d_atoms, natom, atom_types, 
          0.05, rcut, fold, d_mat, d_invmat, k_ewald, n_kvectors, kvectors, norms, cell.volume());

  gpuErrchk(cudaDeviceSynchronize());

  
  for (int i = 0; i < atom_types; i++) {
    for (int j = 0; j < size; j++)
      std::cout << grid[i * size + j] << ' ';
    std::cout << endl;
  }
  cudaFree(grid);
  
  return 0;
}