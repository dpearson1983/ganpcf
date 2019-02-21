#include <iostream>
#include <vector>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector_types.h>
#include "../include/gpuerrchk.h"
#include "../include/ganpcf.hpp"

#define N_threads 64

__constant__ int3 d_shifts[27];
__constant__ float d_R;
__constant__ float d_r;
__constant__ float3 d_L;
__constant__ int d_Nshells;
__constant__ int d_Nparts;

void npcf::initVectors() {
    for (int i = 0; i < npcf::N_shells; ++i) {
        double r1 = npcf::r_min + (i + 0.5)*npcf::Delta_r;
        npcf::rs.push_back(r1);
        npcf::twoPoint.push_back(0.0);
        npcf::DD.push_back(0);
        npcf::DR.push_back(0);
        for (int j = i; j < npcf::N_shells; ++j) {
            double r2 = npcf::r_min + (j + 0.5)*npcf::Delta_r;
            for (int k = j; k < npcf::N_shells; ++k) {
                double r3 = npcf::r_min + (k + 0.5)*npcf::Delta_r;
                if (r3 <= r1 + r2) {
                    npcf::threePoint.push_back(0.0);
                    float3 tri = {(float)r1, (float)r2, (float)r3};
                    npcf::triangles.push_back(tri);
                }
            }
        }
    }
    npcf::DDD.resize(npcf::N_shells*npcf::N_shells*npcf::N_shells);
    npcf::DDR.resize(npcf::N_shells*npcf::N_shells*npcf::N_shells);
    npcf::DRR.resize(npcf::N_shells*npcf::N_shells*npcf::N_shells);
    npcf::RRR.resize(npcf::N_shells*npcf::N_shells*npcf::N_shells);
    
    npcf::w = {0.8888888888888888, 0.5555555555555556, 0.5555555555555556};
    
    npcf::x = {0.0000000000000000, -0.7745966692414834, 0.7745966692414834};
}

void npcf::swapIfGreater(double &a, double &b) {
    if (a > b) {
        double temp = a;
        a = b;
        b = temp;
    }
}

double npcf::sphereOverlapVolume(double d, double R, double r) {
    double V = 0;
    swapIfGreater(r, R);
    if (d < R + r) {
        if (d > R - r) {
            V = (M_PI*(R + r - d)*(R + r - d)*(d*d + 2.0*d*r - 3.0*r*r + 2.0*d*R + 6.0*r*R - 3.0*R*R))/(12.0*d);
        } else {
            V = (4.0*M_PI/3.0)*r*r*r;
        }
    }
    return V;
}

double npcf::crossSectionVolume(double r1, double r2, double r3) {
    double V_oo = sphereOverlapVolume(r1, r3 + 0.5*npcf::Delta_r, r2 + 0.5*npcf::Delta_r);
    double V_oi = sphereOverlapVolume(r1, r3 + 0.5*npcf::Delta_r, r2 - 0.5*npcf::Delta_r);
    double V_io = sphereOverlapVolume(r1, r3 - 0.5*npcf::Delta_r, r2 + 0.5*npcf::Delta_r);
    double V_ii = sphereOverlapVolume(r1, r3 - 0.5*npcf::Delta_r, r2 - 0.5*npcf::Delta_r);
    
    return V_oo - V_oi - V_io + V_ii;
}

int npcf::getPermutations(double r1, double r2, double r3) {
    int perm = 1;
    if (r1 != r2 && r1 != r3 && r2 != r3) {
        perm = 6;
    } else if ((r1 == r2 && r1 != r3) || (r1 == r3 && r1 != r2) || (r2 == r3 && r2 != r1)) {
        perm = 3;
    }
    return perm;
}

double npcf::sphericalShellVolume(double r) {
    double r_o = r + 0.5*npcf::Delta_r;
    double r_i = r - 0.5*npcf::Delta_r;
    return 4.0*M_PI*(r_o*r_o*r_o - r_i*r_i*r_i)/3.0;
}

double npcf::nbarData(double r, double r1) {
    int bin = r/npcf::Delta_r;
    double nbar = npcf::DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
    int num_bins = npcf::DD.size();
    if (r <= (bin + 0.5)*npcf::Delta_r) {
        if (bin != 0) {
            double n1 = npcf::DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = npcf::DD[bin - 1]/(npcf::N_parts*sphericalShellVolume(r1 - npcf::Delta_r));
            double b = n1 - ((n1 - n2)/npcf::Delta_r)*r1;
            nbar = ((n1 - n2)/npcf::Delta_r)*r + b;
        } else {
            double n1 = npcf::DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = npcf::DD[bin + 1]/(npcf::N_parts*sphericalShellVolume(r1 + npcf::Delta_r));
            double b = n1 - ((n2 - n1)/npcf::Delta_r)*r1;
            nbar = ((n2 - n1)/npcf::Delta_r)*r + b;
        }
    } else {
        if (bin != num_bins - 1) {
            double n1 = npcf::DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = npcf::DD[bin + 1]/(npcf::N_parts*sphericalShellVolume(r1 + npcf::Delta_r));
            double b = n1 - ((n2 - n1)/npcf::Delta_r)*r1;
            nbar = ((n2 - n1)/npcf::Delta_r)*r + b;
        } else {
            double n1 = npcf::DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = npcf::DD[bin - 1]/(npcf::N_parts*sphericalShellVolume(r1 - npcf::Delta_r));
            double b = n1 - ((n1 - n2)/npcf::Delta_r)*r1;
            nbar = ((n1 - n2)/npcf::Delta_r)*r + b;
        }
    }
    return nbar;
}

double npcf::gaussQuadCrossSection(double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < npcf::w.size(); ++i) {
        double r_1 = r1 + 0.5*npcf::Delta_r*npcf::x[i];
        result += 0.5*npcf::Delta_r*npcf::w[i]*crossSectionVolume(r_1, r2, r3)*r_1*r_1;
    }
    return result;
}

double npcf::gaussQuadCrossSectionDDR(double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < npcf::w.size(); ++i) {
        double r_1 = r1 + 0.5*npcf::Delta_r*npcf::x[i];
        double nbar = nbarData(r_1, r1);
        result += 0.5*npcf::Delta_r*npcf::w[i]*crossSectionVolume(r_1, r2, r3)*r_1*r_1*nbar;
    }
    return result;
}

npcf::npcf(int timesRandoms, int numShells, double VolBox, double rMin, double rMax) {
    npcf::timesRan = timesRandoms;
    npcf::N_shells = numShells;
    npcf::r_max = rMax;
    npcf::r_min = rMin;
    npcf::V_box = VolBox;
    npcf::Delta_r = (rMax - rMin)/numShells;
    npcf::initVectors();
}

int npcf::getShells(double *shells[]) {
    for (int i = 0; i < npcf::rs.size(); ++i)
        shells[0][i] = npcf::rs[i];
    return 1;
}

int npcf::getNumTriangles() {
    return npcf::triangles.size();
}

int npcf::getTriangles(float3 *tris[]) {
    for (int i = 0; i < npcf::triangles.size(); ++i) {
        tris[0][i].x = npcf::triangles[i].x;
        tris[0][i].y = npcf::triangles[i].y;
        tris[0][i].z = npcf::triangles[i].z;
    }
    return 1;
}
        
int npcf::setNumParticles(int numParticles) {
    npcf::N_parts = numParticles;
    npcf::N_rans = numParticles*npcf::timesRan;
    npcf::nbar_ran = npcf::N_rans/npcf::V_box;
    return 1;
}

void npcf::getRRR() {
    for (int i = 0; i < npcf::N_shells; ++i) {
        for (int j = i; j < npcf::N_shells; ++j) {
            for (int k = j; k < npcf::N_shells; ++k) {
                if (npcf::rs[k] <= npcf::rs[i] + npcf::rs[j]) {
                    int index = k + npcf::N_shells*(j + npcf::N_shells*i);
                    double V = npcf::gaussQuadCrossSection(npcf::rs[i], npcf::rs[j], npcf::rs[k]);
                    int n_perm = npcf::getPermutations(npcf::rs[i], npcf::rs[j], npcf::rs[k]);
                    npcf::RRR[index] = int(4.0*M_PI*n_perm*npcf::nbar_ran*npcf::nbar_ran*V*npcf::N_rans);
                }
            }
        }
    }
}

void npcf::getDRR() {
    for (int i = 0; i < npcf::N_shells; ++i) {
        for (int j = i; j < npcf::N_shells; ++j) {
            for (int k = j; k < npcf::N_shells; ++k) {
                if (npcf::rs[k] <= npcf::rs[i] + npcf::rs[j]) {
                    int index = k + npcf::N_shells*(j + npcf::N_shells*i);
                    double V = npcf::gaussQuadCrossSection(npcf::rs[i], npcf::rs[j], npcf::rs[k]);
                    int n_perm = npcf::getPermutations(npcf::rs[i], npcf::rs[j], npcf::rs[k]);
                    npcf::DRR[index] = int(4.0*M_PI*n_perm*npcf::nbar_ran*npcf::nbar_ran*V*npcf::N_parts);
                }
            }
        }
    }
}

void npcf::getDDR() {
    for (int i = 0; i < npcf::N_shells; ++i) {
        double r1 = npcf::rs[i];
        for (int j = i; j < npcf::N_shells; ++j) {
            double r2 = npcf::rs[j];
            for (int k = j; k < npcf::N_shells; ++k) {
                double r3 = npcf::rs[k];
                if (npcf::rs[k] <= npcf::rs[i] + npcf::rs[j]) {
                   int index = k + npcf::N_shells*(j + npcf::N_shells*i);
                   double V = npcf::gaussQuadCrossSectionDDR(r1, r2, r3);
                   double N_temp = 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                   if (r1 != r2 && r1 != r3 && r2 != r3) {
                       V = npcf::gaussQuadCrossSectionDDR(r2, r3, r1);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(r3, r1, r2);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(r1, r3, r2);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(r2, r1, r3);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(r3, r2, r1);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                   } else if ((r1 == r2 && r1 != r3) || (r1 == r3 && r1 != r2) || (r2 == r3 && r2 != r1)) {
                       V = npcf::gaussQuadCrossSectionDDR(r2, r3, r1);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(r3, r1, r2);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                   }
                   npcf::DDR[index] = int(floor(N_temp + 0.5));
                }
            }
        }
    }
}

void npcf::getDR() {
    for (int i = 0; i < npcf::rs.size(); ++i) {
        npcf::DR[i] = npcf::sphericalShellVolume(npcf::rs[i])*npcf::nbar_ran*npcf::N_parts;
    }
}

int npcf::get2ptSize() {
    return npcf::twoPoint.size();
}

int npcf::get3ptSize() {
    return npcf::threePoint.size();
}

__device__ float get_separation(float3 r1, float3 r2) {
    return sqrtf((r1.x - r2.x)*(r1.x - r2.x) + (r1.y - r2.y)*(r1.y - r2.y) + (r1.z - r2.z)*(r1.z - r2.z));
}

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
    float Delta_r = (d_R - d_r)/d_Nshells;
    int shell1 = int((d1 - d_r)/Delta_r);
    int shell2 = int((d2 - d_r)/Delta_r);
    int shell3 = int((d3 - d_r)/Delta_r);
    return shell3 + d_Nshells*(shell2 + d_Nshells*shell1);
}

__device__ int4 get_index(int4 ngp, int i, int3 n, float3 &rShift) {
    ngp.x += d_shifts[i].x;
    ngp.y += d_shifts[i].y;
    ngp.z += d_shifts[i].z;
    rShift.x = 0.0;
    rShift.y = 0.0;
    rShift.z = 0.0;
    if (ngp.x == n.x) {
        ngp.x = 0;
        rShift.x = d_L.x;
    }
    if (ngp.y == n.y) {
        ngp.y = 0;
        rShift.y = d_L.y;
    }
    if (ngp.z == n.z) {
        ngp.z = 0;
        rShift.z = d_L.z;
    }
    if (ngp.x == -1) {
        ngp.x = n.x - 1;
        rShift.x = -d_L.x;
    }
    if (ngp.y == -1) {
        ngp.y = n.y - 1;
        rShift.y = -d_L.y;
    }
    if (ngp.z == -1) {
        ngp.z = n.z - 1;
        rShift.z = -d_L.z;
    }
    ngp.w = ngp.z + n.z*(ngp.y + n.y*ngp.x);
    return ngp;
}

__global__ void countPairs(float3 *d_p1, float3 **d_p2, int *p2_sizes, int *d_partsPerShell, int3 n) {
    // Calculate the thread ID for the current GPU thread
    int tid = threadIdx.x + blockIdx.x*blockDim.x;
    float Delta_r = (d_R - d_r)/d_Nshells;
    
    if (tid < d_Nparts) {
        float3 r1 = d_p1[tid];
        int4 ngp1 = {int(r1.x/d_R), int(r1.y/d_R), int(r1.z/d_R), 0};
        for (int i = 0; i < 27; ++i) {
            float3 rShift2;
            int4 index = get_index(ngp1, i, n, rShift2);
            int size2 = p2_sizes[index.w];
            for (int part2 = 0; part2 < size2; ++part2) {
                float3 r2 = d_p2[index.w][part2];
                r2.x += rShift2.x;
                r2.y += rShift2.y;
                r2.z += rShift2.z;
                float dist = get_separation(r1, r2);
                if (dist < d_R && dist > d_r) {
                    int shell = int((dist - d_r)/Delta_r);
                    atomicAdd(&d_partsPerShell[shell], 1);
                }
            }
        }
    }
}

__global__ void countTriangles(float3 *d_p1, float3 **d_p2, float3 **d_p3, int *p2_sizes, int *p3_sizes, 
                               int *d_triangles, int3 n) {
    // Calculate the thread ID for the current GPU thread
    int tid = threadIdx.x + blockIdx.x*blockDim.x;
    
    if (tid < d_Nparts) {
        float3 r1 = d_p1[tid];
        int4 ngp1 = {int(r1.x/d_R), int(r1.y/d_R), int(r1.z/d_R), 0};
        if (ngp1.x == n.x) ngp1.x--;
        if (ngp1.y == n.x) ngp1.y--;
        if (ngp1.z == n.x) ngp1.z--;
        for (int i = 0; i < 27; ++i) {
            float3 rShift2;
            int4 index = get_index(ngp1, i, n, rShift2);
            int size2 = p2_sizes[index.w];
            for (int part2 = 0; part2 < size2; ++part2) {
                float3 r2 = d_p2[index.w][part2];
                r2.x += rShift2.x;
                r2.y += rShift2.y;
                r2.z += rShift2.z;
                float d1 = get_separation(r1, r2);
                if (d1 < d_R && d1 > d_r) {
                    for (int j = 0; j < 27; ++j) {
                        float3 rShift3;
                        int4 index2 = get_index(ngp1, j, n, rShift3);
                        int size3 = p3_sizes[index2.w];
                        for (int part3 = 0; part3 < size3; ++part3) {
                            float3 r3 = d_p3[index2.w][part3];
                            r3.x += rShift3.x;
                            r3.y += rShift3.y;
                            r3.z += rShift3.z;
                            float d2 = get_separation(r1, r3);
                            float d3 = get_separation(r2, r3);
                            if (d2 < d_R && d3 < d_R && d2 > d_r && d3 > d_r) {
                                int shell = get_shell(d1, d2, d3);
                                atomicAdd(&d_triangles[shell], 1);
                            }
                        }
                    }
                }
            }
        }
    }
}

std::vector<int3> getShifts() {
    std::vector<int3> shifts;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                int3 temp = {i, j, k};
                shifts.push_back(temp);
            }
        }
    }
    return shifts;
}

void npcf::rezeroVectors() {
#pragma omp parallel for
    for (int i = 0; i < npcf::DD.size(); ++i) {
        npcf::DD[i] = 0;
        npcf::DR[i] = 0;
    }
#pragma omp parallel for
    for (int i = 0;  i < npcf::DDD.size(); ++i) {
        npcf::DDD[i] = 0;
        npcf::DDR[i] = 0;
        npcf::DRR[i] = 0;
        npcf::RRR[i] = 0;
    }
}

int npcf::calculateCorrelations(float3 *galaxies[]) {
    float l = (float)pow(npcf::V_box, 1.0/3.0);
    float3 L = {l, l, l};
    float R = (float)npcf::r_max;
    float r = (float)npcf::r_min;
    std::vector<int3> shifts = getShifts();
    
    npcf::rezeroVectors();
    
    int numParts = npcf::N_parts;
    int numShells = npcf::N_shells;
    
    gpuErrchk(cudaMemcpyToSymbol(d_shifts, shifts.data(), shifts.size()*sizeof(int3)));
    gpuErrchk(cudaMemcpyToSymbol(d_L, &L, sizeof(float3)));
    gpuErrchk(cudaMemcpyToSymbol(d_R, &R, sizeof(float)));
    gpuErrchk(cudaMemcpyToSymbol(d_r, &r, sizeof(float)));
    gpuErrchk(cudaMemcpyToSymbol(d_Nparts, &numParts, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(d_Nshells, &numShells, sizeof(int)));
    
    // Find the cell size such that the box size is an integer multiple of the cell size
    int l_c = floor(L.x/R);
    double L_c = L.x/l_c;
    int3 N = {l_c, l_c, l_c};
    
    std::vector<std::vector<float3>> gals(N.x*N.y*N.z);
    std::vector<int> sizes;
    float3 **d_gals;
    float3 *d_galaxies;
    int *d_DD, *d_DDD, *d_sizes;
    
    for (int i = 0; i < npcf::N_parts; ++i) {
        int ix = galaxies[0][i].x/L_c;
        if (ix == N.x) ix--;
        int iy = galaxies[0][i].y/L_c;
        if (iy == N.y) iy--;
        int iz = galaxies[0][i].z/L_c;
        if (iz == N.z) iz--;
        int index = iz + N.z*(iy + N.y*ix);
        gals[index].push_back(galaxies[0][i]);
    }
    
    std::vector<float3> gal_vec;
    float3 **h_gals = (float3 **)malloc(gals.size()*sizeof(float3 *));
    for (int i = 0; i < gals.size(); ++i) {
        sizes.push_back(gals[i].size());
        gpuErrchk(cudaMalloc((void **)&h_gals[i], gals[i].size()*sizeof(float3)));
        gpuErrchk(cudaMemcpy(h_gals[i], gals[i].data(), gals[i].size()*sizeof(float3), cudaMemcpyHostToDevice));
        for (int j = 0; j < gals[i].size(); ++j) {
//             galaxies[0][j].x = gals[i][j].x;
//             galaxies[0][j].y = gals[i][j].y;
//             galaxies[0][j].z = gals[i][j].z;
            gal_vec.push_back(gals[i][j]);
        }
    }
    gpuErrchk(cudaMalloc(&d_gals, gals.size()*sizeof(float3 *)));
    gpuErrchk(cudaMemcpy(d_gals, h_gals, gals.size()*sizeof(float3 *), cudaMemcpyHostToDevice));
    
    gpuErrchk(cudaMalloc((void **)&d_galaxies, npcf::N_parts*sizeof(float3)));
    gpuErrchk(cudaMalloc((void **)&d_sizes, sizes.size()*sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&d_DD, npcf::DD.size()*sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&d_DDD, npcf::DDD.size()*sizeof(int)));
    
    gpuErrchk(cudaMemcpy(d_galaxies, gal_vec.data(), gal_vec.size()*sizeof(float3), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_sizes, sizes.data(), sizes.size()*sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_DD, npcf::DD.data(), npcf::DD.size()*sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_DDD, npcf::DDD.data(), npcf::DDD.size()*sizeof(int), cudaMemcpyHostToDevice));
    
    int N_blocks = npcf::N_parts/N_threads + 1;
    
    countPairs<<<N_blocks, N_threads>>>(d_galaxies, d_gals, d_sizes, d_DD, N);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    countTriangles<<<N_blocks, N_threads>>>(d_galaxies, d_gals, d_gals, d_sizes, d_sizes, d_DDD, N);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    
    gpuErrchk(cudaMemcpy(npcf::DD.data(), d_DD, npcf::DD.size()*sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(npcf::DDD.data(), d_DDD, npcf::DDD.size()*sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaDeviceSynchronize());
    
    getDR();
    getRRR();
    getDRR();
    getDDR();
    
    double alpha = double(npcf::N_parts)/double(npcf::N_rans);
    std::cout << alpha << std::endl;
    for (int i = 0; i < npcf::twoPoint.size(); ++i) {
        std::cout << npcf::DD[i] << std::endl;
        npcf::twoPoint[i] = double(npcf::DD[i])/(alpha*double(npcf::DR[i])) - 1.0;
    }
    
    for (int i = 0; i < npcf::threePoint.size(); ++i) {
        int ix = npcf::triangles[i].x/npcf::Delta_r;
        int iy = npcf::triangles[i].y/npcf::Delta_r;
        int iz = npcf::triangles[i].z/npcf::Delta_r;
        int index = iz + npcf::N_shells*(iy + npcf::N_shells*ix);
        npcf::threePoint[i] = (double(npcf::DDD[index]) - 3.0*alpha*double(npcf::DDR[index]) 
                               + 3.0*alpha*alpha*double(npcf::DRR[index]) -
                               alpha*alpha*alpha*double(npcf::RRR[index]))/
                               (alpha*alpha*alpha*double(npcf::RRR[index]));
    }
    
    cudaFree(d_DDD);
    cudaFree(d_DD);
    cudaFree(d_sizes);
    cudaFree(d_galaxies);
    for (int i = 0; i < gals.size(); ++i) {
        cudaFree(h_gals[i]);
    }
    cudaFree(d_gals);
    delete[] h_gals;
    return 1;
}

int npcf::get2pt(double *twoPt[]) {
    for (int i = 0; i < npcf::twoPoint.size(); ++i) {
        twoPt[0][i] = npcf::twoPoint[i];
    }
    return 1;
}

int npcf::get3pt(double *threePt[]) {
    for (int i = 0; i < npcf::threePoint.size(); ++i) {
        threePt[0][i] = npcf::threePoint[i];
    }
    return 1;
}
