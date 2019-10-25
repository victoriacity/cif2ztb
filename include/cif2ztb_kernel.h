#ifndef CIFZTBH
#define CIFZTBH

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "vec3_cuda.h"
#include "unitcell.h"
using namespace std;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__host__ int get_kvectors(vec3*, float*, float, float*, float*);

__global__
void create_grid_pos(vec3*, vec3, int[], bool, vec3);

__global__
void calculate_grid(vec3*, uint32_t, vec3*, int*, int, int, float, float, bool, 
                float*, float*, float, int, vec3*, float*, float);

const uint32_t MAX_KVECTORS = 10000;
const int MCCCS_EPS = 1e6;
const float RCUTMIN = 0.05;
const bool FOLD = true;
const uint32_t CUDA_BLOCKSIZE = 256;

#endif