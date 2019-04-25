#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <cmath>
#include <cuda.h>
#include <omp.h>
#include <vector_types.h>
#include <harppi.h>
#include "gpuerrchk.h"

// Set some values for easy reference and updating
#define PI 3.1415926535897932384626433832795
#define N_threads 256

__constant__ float d_L;
__constant__ float d_R;
__constant__ int d_Nshells;
__constant__ int d_Nparts;

// Declare arrays or variables that can be stored in the GPU's constant memory. These cannot be changed during
// the execution of the GPU kernel, but can be updated repeated from the host side code.
__constant__ int3 d_shifts[27];

// The __device__ decorator specifies a function that can be called from inside the __global__ decorated GPU
// kernel, or any other __device__ decorated function. The decorators lets nvcc know which parts of the code
// it needs to compile, while all other code is passed along to a standard compiler like the  GNU Compiler
// collection's (GCC) g++. As the name of this function suggests, it simply calculates the 3D separation of two
// points.
__device__ float get_separation(float3 &r1, float3 &r2) {
    return sqrtf((r1.x - r2.x)*(r1.x - r2.x) + (r1.y - r2.y)*(r1.y - r2.y) + (r1.z - r2.z)*(r1.z - r2.z));
}

// The code makes sure that the three lengths can in fact form a triangle. The __host__ decorator means that this
// can also be called from the CPU.
__device__ __host__ bool isTriangle(float r1, float r2, float r3) {
    if (r1 > r2) {
        float temp = r1;
        r1 = r2;
        r2 = temp;
    }
    if (r1 > r3) {
        float temp = r1;
        r1 = r3;
        r3 = temp;
    }
    if (r2 > r3) {
        float temp = r2;
        r2 = r3;
        r3 = temp;
    }
    if (r3 <= r1 + r2) {
        return true;
    } else {
        return false;
    }
}

// Determines which bin a triangle belongs in. By sorting the lengths first, we can make sure that all permutations
// of lengths are put in the same bin.
__device__ int get_shell(float d1, float d2, float d3) {
    if (d1 > d2) {
        float temp = d1;
        d1 = d2;
        d2 = temp;
    }
    if (d1 > d3) {
        float temp = d1;
        d1 = d3;
        d3 = temp;
    }
    if (d2 > d3) {
        float temp = d2;
        d2 = d3;
        d3 = temp;
    }
    if (d1 <= d2 && d2 <= d3 && d3 <= d1 + d2) {
        int shell1 = d1*d_Nshells/d_R;
        int shell2 = d2*d_Nshells/d_R;
        int shell3 = d3*d_Nshells/d_R;
        return shell3 + d_Nshells*(shell2 + d_Nshells*shell1);
    } else {
        return 1.0;
    }
}

// The points are binned so that points that are spatially close are stored near each other. This function
// is used to cycle through all the bins that directly neighbor the bin containing the first point, including the 
// bin that holds that first point, to find triangle. This way, even though the algorithm is still O(N^3), the N
// can be much smaller.
__device__ int4 get_index(int4 ngp, int i, int n, float3 &rShift) {
    ngp.x += d_shifts[i].x;
    ngp.y += d_shifts[i].y;
    ngp.z += d_shifts[i].z;
    rShift.x = 0.0;
    rShift.y = 0.0;
    rShift.z = 0.0;
    if (ngp.x >= n) {
        ngp.x -= n;
        rShift.x = d_L;
    }
    if (ngp.y >= n) {
        ngp.y -= n;
        rShift.y = d_L;
    }
    if (ngp.z >= n) {
        ngp.z -= n;
        rShift.z = d_L;
    }
    if (ngp.x <= -1) {
        ngp.x += n;
        rShift.x = -d_L;
    }
    if (ngp.y <= -1) {
        ngp.y += n;
        rShift.y = -d_L;
    }
    if (ngp.z <= -1) {
        ngp.z += n;
        rShift.z = -d_L;
    }
    ngp.w = ngp.z + n*(ngp.y + n*ngp.x);
    return ngp;
}

__global__ void countPairs(float3 *d_p, int2 *d_cells, int *d_pairs, int n) {
    int tid = threadIdx.x + blockIdx.x*blockDim.x;
    
    if (tid < d_Nparts) {
        float3 p1 = d_p[tid];
        int4 ngp1 = {int(p1.x/d_R), int(p1.y/d_R), int(p1.z/d_R), 0};
        if (ngp1.x == n) ngp1.x--;
        if (ngp1.y == n) ngp1.y--;
        if (ngp1.z == n) ngp1.z--;
        for (int i = 0; i < 27; ++i) {
            float3 p2shift;
            int4 index2 = get_index(ngp1, i, n, p2shift);
            int2 bounds = d_cells[index2.w];
            for (int part2 = bounds.x; part2 <= bounds.y; ++part2) {
                float3 p2 = d_p[part2];
                p2.x += p2shift.x;
                p2.y += p2shift.y;
                p2.z += p2shift.z;
                float r1 = get_separation(p1, p2);
                if (r1 < d_R && r1 > 0) {
                    int shell = int(r1*d_Nshells/d_R);
                    atomicAdd(&d_pairs[shell], 1);
                }
            }
        }
    }
}

__global__ void countTriangles(float3 *d_p, int2 *d_cells, int *d_triangles, int n) {
    int tid = threadIdx.x + blockIdx.x*blockDim.x;
    
    if (tid < d_Nparts) {
        float3 p1 = d_p[tid];
        int4 ngp1 = {int(p1.x/d_R), int(p1.y/d_R), int(p1.z/d_R)};
        if (ngp1.x == n) ngp1.x--;
        if (ngp1.y == n) ngp1.y--;
        if (ngp1.z == n) ngp1.z--;
        for (int i = 0; i < 27; ++i) {
            float3 p2shift;
            int4 index2 = get_index(ngp1, i, n, p2shift);
            int2 bounds2 = d_cells[index2.w];
            for (int part2 = bounds2.x; part2 <= bounds2.y; ++part2) {
                float3 p2 = d_p[part2];
                p2.x += p2shift.x;
                p2.y += p2shift.y;
                p2.z += p2shift.z;
                float r1 = get_separation(p1, p2);
                if (r1 < d_R && r1 > 0) {
                    for (int j = 0; j < 27; ++j) {
                        float3 p3shift;
                        int4 index3 = get_index(ngp1, j, n, p3shift);
                        int2 bounds3 = d_cells[index3.w];
                        for (int part3 = bounds3.x; part3 <= bounds3.y; ++part3) {
                            float3 p3 = d_p[part3];
                            p3.x += p3shift.x;
                            p3.y += p3shift.y;
                            p3.z += p3shift.z;
                            float r2 = get_separation(p2, p3);
                            float r3 = get_separation(p1, p3);
                            if (r2 < d_R && r2 > 0 && r3 < d_R && r3 > 0) {
                                int shell = get_shell(r1, r2, r3);
                                atomicAdd(&d_triangles[shell], 1);
                            }
                        }
                    }
                }
            }
        }
    }
}

std::vector<int2> getCells(std::vector<float3> &parts, double L, double R, int &n) {
    n = int(L/R);
    std::vector<std::vector<float3>> H(n*n*n);
    std::vector<int2> cells;
    
    std::cout << "Binning particles..." << std::endl;
    for (int i = 0; i < parts.size(); ++i) {
        int ix = parts[i].x/R;
        int iy = parts[i].y/R;
        int iz = parts[i].z/R;
        int index = iz + n*(iy + n*ix);
        H[index].push_back(parts[i]);
    }
    
    int part = 0;
    for (int i = 0; i < H.size(); ++i) {
        int2 cell = {part, int(part + H[i].size() - 1)};
        cells.push_back(cell);
        for (int j = 0; j < H[i].size(); ++j) {
            parts[part + j] = H[i][j];
        }
        part += H[i].size();
    }
    
    return cells;
}

void writePairs(std::string file, std::vector<int> &pairs, float R, int N_shells) {
    double dr = R/N_shells;
    std::ofstream fout(file);
    for (int i = 0; i < pairs.size(); ++i) {
        double r = (i + 0.5)*dr;
        fout << r << " " << pairs[i] << "\n";
    }
    fout.close();
}

void writeTriangles(std::string file, std::vector<int> &triangles, float R, int N_shells) {
    double dr = R/N_shells;
    std::ofstream fout(file);
    for (int i = 0; i < N_shells; ++i) {
        double r1 = (i + 0.5)*dr;
        for (int j = i; j < N_shells; ++j) {
            double r2 = (j + 0.5)*dr;
            for (int k = j; k < N_shells; ++k) {
                double r3 = (k + 0.5)*dr;
                if (isTriangle(r1, r2, r3)) {
                    int index = k + N_shells*(j + N_shells*i);
                    fout << r1 << " " << r2 << " " << r3 << " " << triangles[index] << "\n";
                }
            }
        }
    }
    fout.close();
}

int main(int argc, char *argv[]) {
    parameters p(argv[1]);
    p.print();
    
    std::vector<int3> shifts;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                int3 temp = {i, j, k};
                shifts.push_back(temp);
//                 std::cout << i << ", " << j << ", " << k << std::endl;
            }
        }
    }
    
    float L = float(p.getd("L"));
    float R = float(p.getd("R"));
    int N_shells = p.geti("N_shells");
    
    cudaSetDevice(0);
    
    std::cout << "Writing values to constant memory..." << std::endl;
    gpuErrchk(cudaMemcpyToSymbol(d_L, &L, sizeof(float)));
    gpuErrchk(cudaMemcpyToSymbol(d_R, &R, sizeof(float)));
    gpuErrchk(cudaMemcpyToSymbol(d_Nshells, &N_shells, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(d_shifts, shifts.data(), shifts.size()*sizeof(int3)));
    
    std::cout << "Reading input file..." << std::endl;
    int N_parts;
    std::ifstream fin(p.gets("inFile"), std::ios::binary);
    fin.read((char *)&N_parts, sizeof(int));
    std::cout << N_parts << std::endl;
    std::vector<float3> parts(N_parts);
    fin.read((char *)parts.data(), parts.size()*sizeof(float3));
    fin.close();
    
    std::cout << "Writing number of particles to constant memory..." << std::endl;
    std::cout << "    Number or particles = " << N_parts << std::endl;
    cudaMemcpyToSymbol(d_Nparts, &N_parts, sizeof(int));
    
    std::cout << "Setting up cells..." << std::endl;
    int n;
    std::vector<int2> cells = getCells(parts, p.getd("L"), p.getd("R"), n);
    std::ofstream fout("cells.dat");
    for (int i = 0; i < cells.size(); ++i) {
        fout << cells[i].x << " " << cells[i].y << "\n";
    }
    fout.close();
    std::vector<int> pairs(N_shells);
    std::vector<int> triangles(N_shells*N_shells*N_shells);
    
    std::cout << "Declaring device pointers..." << std::endl;
    int2 *d_cells;
    int *d_pairs, *d_triangles;
    float3 *d_parts;
    
    std::cout << "Allocating device pointers..." << std::endl;
    gpuErrchk(cudaMalloc((void **)&d_cells, cells.size()*sizeof(int2)));
    gpuErrchk(cudaMalloc((void **)&d_pairs, pairs.size()*sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&d_triangles, triangles.size()*sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&d_parts, parts.size()*sizeof(float3)));
    
    std::cout << "Copying from device to host..." << std::endl;
    gpuErrchk(cudaMemcpy(d_cells, cells.data(), cells.size()*sizeof(int2), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_pairs, pairs.data(), pairs.size()*sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_triangles, triangles.data(), triangles.size()*sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_parts, parts.data(), parts.size()*sizeof(float3), cudaMemcpyHostToDevice));
    
    std::cout << "Executing GPU kernels..." << std::endl;
    int N_blocks = N_parts/N_threads + 1;
    cudaEvent_t begin, end;
    float elapsedTime;
    cudaEventCreate(&begin);
    cudaEventRecord(begin, 0);
    countPairs<<<N_blocks, N_threads>>>(d_parts, d_cells, d_pairs, n);
    cudaEventCreate(&end);
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    cudaEventElapsedTime(&elapsedTime, begin, end);
    std::cout << "Time to count pairs: " << elapsedTime << " ms" << std::endl;
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    cudaEventCreate(&begin);
    cudaEventRecord(begin, 0);
    countTriangles<<<N_blocks, N_threads>>>(d_parts, d_cells, d_triangles, n);
    cudaEventCreate(&end);
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    cudaEventElapsedTime(&elapsedTime, begin, end);
    std::cout << "Time to count triangles: " << elapsedTime << " ms" << std::endl;
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    
    std::cout << "Reading data back from GPU..." << std::endl;
    gpuErrchk(cudaMemcpy(pairs.data(), d_pairs, pairs.size()*sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(triangles.data(), d_triangles, triangles.size()*sizeof(int), cudaMemcpyDeviceToHost));
    
    std::cout << "Writing output files..." << std::endl;
    writePairs(p.gets("pairsFile"), pairs, R, N_shells);
    writeTriangles(p.gets("triangleFile"), triangles, R, N_shells);
    
    std::cout << "Freeing GPU memory..." << std::endl;
    gpuErrchk(cudaFree(d_cells));
    gpuErrchk(cudaFree(d_pairs));
    gpuErrchk(cudaFree(d_triangles));
    gpuErrchk(cudaFree(d_parts));
    
    std::cout << "Done!" << std::endl;
    return 0;
}
