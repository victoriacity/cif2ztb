#include "cif2ztb_kernel.h"


void cif2ztb(string in_file, string out_file, 
             float rcut, float k_ewald, float spacing, bool use_fractional_basis) {
  
  UnitCell cell_min, cell;
  cell_min.read_from_cif(in_file);
  
  int* dims;
  vec3 duplicates;
  gpuErrchk(cudaMallocManaged(&dims, 3 * sizeof(int)));
  for (int i = 0; i < 3; i++) {
    dims[i] = (int) round(cell_min.cell_length[i] / spacing);
    duplicates[i] = (int) ceil(2 * rcut / cell_min.cell_length[i]);
  }
  cell = cell_min.duplicate(duplicates);
  
  int atom_types;
  uint32_t natom, size;  
  vector<vec3> coordinates;
  vector<int> atoms;
  cell.get_all_coordinates(&coordinates, &atoms, false);

  atom_types = cell.atomtypes.size();
  natom = coordinates.size();
  size = dims[0] * dims[1] * dims[2];
  
 
  float* d_mat, *invmat, *d_invmat, *norms;
  vec3* d_coordinates, *kvectors;
  int* d_atoms;
  int n_kvectors;

  gpuErrchk(cudaMalloc(&d_mat, 9 * sizeof(float)));
  gpuErrchk(cudaMalloc(&d_invmat, 9 * sizeof(float)));
  invmat = invert_upper_trig_matrix(cell.transformation_matrix());
  gpuErrchk(cudaMemcpy(d_mat, cell.transformation_matrix(), 9 * sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_invmat, invmat, 9 * sizeof(float), cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(&d_coordinates, natom * sizeof(vec3)));
  gpuErrchk(cudaMalloc(&d_atoms, natom * sizeof(int)));
  gpuErrchk(cudaMemcpy(d_coordinates, coordinates.data(), natom * sizeof(vec3), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_atoms, atoms.data(), natom * sizeof(int), cudaMemcpyHostToDevice));

  if (k_ewald > 0) {
    gpuErrchk(cudaMallocManaged(&norms, MAX_KVECTORS * sizeof(float)));
    gpuErrchk(cudaMallocManaged(&kvectors, MAX_KVECTORS * sizeof(vec3)));
    n_kvectors = get_kvectors(kvectors, norms, k_ewald, cell.transformation_matrix(), invmat);
  }

  vec3* grid;
  gpuErrchk(cudaMallocManaged(&grid, atom_types * size * sizeof(vec3)));

  uint32_t blockSize = CUDA_BLOCKSIZE;
  uint32_t numBlocks = (size + blockSize - 1) / blockSize;

  create_grid_pos<<<numBlocks, blockSize>>>(grid, spacing, dims, use_fractional_basis, cell.cell_angle);

  calculate_grid<<<numBlocks, blockSize>>>(grid, size, d_coordinates, d_atoms, natom, atom_types, 
          RCUTMIN, rcut, FOLD, d_mat, d_invmat, k_ewald, n_kvectors, kvectors, norms, cell.volume());

  gpuErrchk(cudaDeviceSynchronize());

  ofstream file;
  file.open(out_file);
  
  file << cell_min.cell_length << ' ' << cell_min.cell_angle << ' ' << spacing << endl;
  file << rcut << ' ' << k_ewald << ' ' << (int) use_fractional_basis << endl;

  for (int i = 0; i < atom_types; i++) {
    for (int j = 0; j < size; j++)
      file << grid[i * size + j] << ' ';
    file << endl;
  }

  cudaFree(grid);
  file.close();
}