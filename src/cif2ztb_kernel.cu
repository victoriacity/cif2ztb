#include "cif2ztb_kernel.h"


void cif2ztb(string in_file, string out_file, 
             float rcut, float k_ewald, float spacing, bool use_fractional_basis, bool out_binary) {
  // overrides fractional basis
  if (out_binary) use_fractional_basis = true;
  UnitCell cell_min, cell;
  cell_min.read_from_cif(in_file);
  int* dims;
  vec3 duplicates, spacing_all;
  vec3 dim_mask[3] = {vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1)};
  gpuErrchk(cudaMallocManaged(&dims, 3 * sizeof(int)));
  for (int i = 0; i < 3; i++) {
    // Also use soft spacing if using fractional basis
    // to be consistent with MCCCS-MN
    if (use_fractional_basis) {
      dims[i] = (int) (cell_min.cell_length[i] / spacing);
      spacing_all[i] = cell_min.cell_length[i] / dims[i];
    }
    else {
      dims[i] = (int) round(
        cell_min.frac_to_abs(cell_min.cell_length * dim_mask[i])[i] 
          / spacing);
      spacing_all[i] = spacing;
    }
    duplicates[i] = (int) ceil(2 * rcut / cell_min.frac_to_abs(cell_min.cell_length * dim_mask[i])[i]);
  }
  vec3 duplicates_tmp = duplicates;
  // reads in the first line of CIF to get number of unit cells
  ifstream file_cif(in_file);
  if (file_cif.is_open()) {
    string firstline;
    getline(file_cif, firstline);
    if (firstline[0] == '#') {
      firstline = firstline.substr(1);
      stringstream ss;
      ss << firstline;
      for (int i = 0; i < 3; i++) {
        ss >> duplicates[i];
      }
    }
    file_cif.close();
  }

  /*
  for (int i = 0; i < 3; i++) {
    if (duplicates[i] != duplicates_tmp[i]) {
      cout << "Warning: inconsitent cell numbers in " << in_file << ", ";
      cout << "calculated: " << duplicates_tmp << ", read: " << duplicates << endl;
      exit(1);
    }
  }
  exit(0);
  */
  cell = cell_min.duplicate(duplicates);

  int atom_types;
  uint32_t natom, size;  
  vector<vec3> coordinates;
  vector<int> atoms;
  cell.get_all_coordinates(&coordinates, &atoms, false);

  atom_types = cell.atomtypes.size();
  natom = coordinates.size();
  size = dims[0] * dims[1] * dims[2];

  // calculate chemical forula
  int atom_count[atom_types] = {0};
  int dup_total = (int) (duplicates[0] * duplicates[1] * duplicates[2]);
  for (int i = 0; i < natom; i++) {
    atom_count[atoms[i]]++;
  }
  cout << "Chemical formula: ";
  for (int i = 0; i < atom_types; i++) {
    cout << cell.atomtypes[i] << atom_count[i] / dup_total;
  }
  cout << endl << "Duplicates:" << duplicates << endl;

  //for (vec3 r : cell_min.coordinates) cout << r << endl;

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

  create_grid_pos<<<numBlocks, blockSize>>>(grid, spacing_all, dims, use_fractional_basis, cell.cell_angle);

  calculate_grid<<<numBlocks, blockSize>>>(grid, size, d_coordinates, d_atoms, natom, atom_types, 
          RCUTMIN, rcut, FOLD, d_mat, d_invmat, k_ewald, n_kvectors, kvectors, norms, cell.volume());

  gpuErrchk(cudaDeviceSynchronize());

  ofstream file;
  if (out_binary) {
    file.open(out_file, ios::out | ios::binary);
    // Cell parameters, 3*FLOAT64
    for (int i = 0; i < 3; i++) {
      double val = (double) cell_min.cell_length[i];
      val = round(val * MCCCS_EPS ) / MCCCS_EPS;
      file.write((char*) &val, sizeof(double));
    }
    // Cell angles, 3*FLOAT64
    for (int i = 0; i < 3; i++) {
      double val = (double) cell_min.cell_angle[i] * PI / 180;
      file.write((char*) &val, sizeof(double));
    }
    // grid dimensions, 3*INT32
    for (int i = 0; i < 3; i++) {
      file.write((char*) &dims[i], sizeof(int32_t));
    }
    // number of atom types, 1*INT32
    file.write((char*) &atom_types, sizeof(int32_t));
    // Use ewald summation, 1*INT32
    int32_t use_ewald = k_ewald > 0 ? 1 : 0;
    file.write((char*) &use_ewald, sizeof(int32_t));
    // Use tail correction, 1*INT32, always FALSE
    int32_t _false = 0;
    file.write((char*) &_false, sizeof(int32_t));
    // Use shifted potential, 1*INT32, always FALSE
    file.write((char*) &_false, sizeof(int32_t));
    // cutoff radius, 1*FLOAT64
    double rcut_d = (double) rcut;
    file.write((char*) &(rcut_d), sizeof(double));
    // Atom names, 128*CHAR
    for (string atom : cell.atomtypes) {
      char chars[128] = {' '};
      strcpy(chars, atom.c_str());
      for (int i = 0; i < 128; i++) {
        if (chars[i] == '\x00') chars[i] = ' ';
      }
      file.write((char*) chars, 128 * sizeof(char));
    }
    // Write tabulated potential, use FORTRAN order
    // Loop order: atom - channel - k - j - i
    for (int k = 0; k < dims[2]; k++) {
      for (int j = 0; j < dims[1]; j++) {
        for (int i = 0; i < dims[0]; i++) {
          for (int i_atom = 0; i_atom < atom_types; i_atom++) {
            for (int d = 0; d < 3; d++) {
              uint32_t index_lin = i * dims[1] * dims[2] + j * dims[2] + k;
              double val = (double) grid[i_atom * size + index_lin][d];
              file.write((char*) &val, sizeof(double));
            }
          }
        }
      }
    }


  } else {
    file.open(out_file);
    file << dims[0] << ' ' << dims[1] << ' ' << dims[2];
    file << ' ' << cell_min.cell_angle << ' ' << endl;
    file << spacing << ' ' << rcut << ' ' << k_ewald << ' ' << (int) use_fractional_basis << endl;
    for (string atom : cell.atomtypes) {
      file << atom << ' ';
    }
    file << endl;

    for (int i = 0; i < atom_types; i++) {
      for (uint32_t j = 0; j < size; j++)
        file << grid[i * size + j] << ' ';
      file << endl;
    }
  }
  cudaFree(grid);
  file.close();
}