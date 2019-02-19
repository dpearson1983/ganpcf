#include <vector>
#include <cmath>
#include <cuda.h>
#include <vector_types.h>
#include "../include/ganpcf.hpp"

__constant__ int3 d_shifts[27];
__constant__ float d_R;
__constant__ float d_r;
__constant__ float d_L;
__constant__ d_Nshells;
__constant__ d_Nparts;

void npcf::initVectors() {
    for (int i = 0; i < npcf::N_shells; ++i) {
        double r1 = npcf::r_min + (i + 0.5)*npcf::Delta_r;
        npcf::rs.push_back(r1);
        npcf::twoPoint.push_back(0.0);
        npcf::DD.push_back(0);
        npcf::DR.push_back(0);
        for (int j = i; j < npcf::N_shells; ++j) {
            double r2 = npcf::r_min + (j + 0.5)*npcf::Delta_r
            for (int k = j; k < npcf::N_shells; ++k) {
                double r3 = npcf::r_min + (k + 0.5)*npcf::Delta_r
                if (r3 <= r1 + r2) {
                    npcf::threePoint.push_back(0.0);
                    npcf::DDD.push_back(0);
                    npcf::DDR.push_back(0);
                    npcf::DRR.push_back(0);
                    npcf::RRR.push_back(0);
                    float3 tri = {(float)r1, (float)r2, (float)r3};
                    npcf::triangles.push_back(tri);
                }
            }
        }
    }
    
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

double npcf::nbarData(std::vector<int> &DD, double r, double r1) {
    int bin = r/npcf::Delta_r;
    double nbar = DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
    int num_bins = DD.size();
    if (r <= (bin + 0.5)*npcf::Delta_r) {
        if (bin != 0) {
            double n1 = DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin - 1]/(npcf::N_parts*sphericalShellVolume(r1 - npcf::Delta_r));
            double b = n1 - ((n1 - n2)/npcf::Delta_r)*r1;
            nbar = ((n1 - n2)/npcf::Delta_r)*r + b;
        } else {
            double n1 = DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin + 1]/(npcf::N_parts*sphericalShellVolume(r1 + npcf::Delta_r));
            double b = n1 - ((n2 - n1)/npcf::Delta_r)*r1;
            nbar = ((n2 - n1)/npcf::Delta_r)*r + b;
        }
    } else {
        if (bin != num_bins - 1) {
            double n1 = DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin + 1]/(npcf::N_parts*sphericalShellVolume(r1 + npcf::Delta_r));
            double b = n1 - ((n2 - n1)/npcf::Delta_r)*r1;
            nbar = ((n2 - n1)/npcf::Delta_r)*r + b;
        } else {
            double n1 = DD[bin]/(npcf::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin - 1]/(npcf::N_parts*sphericalShellVolume(r1 - npcf::Delta_r));
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

double npcf::gaussQuadCrossSectionDDR(std::vector<int> &DD, double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < npcf::w.size(); ++i) {
        double r_1 = r1 + 0.5*npcf::Delta_r*npcf::x[i];
        double nbar = nbarData(DD, r_1, r1);
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
    npcf::nbar_ran = numRandoms/VolBox;
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
        
int npcf::setNumParts(int numParticles) {
    npcf::N_parts = numParticles;
    npcf::N_rans = numParticles*npcf::timesRan;
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
                   double V = npcf::gaussQuadCrossSectionDDR(DD, r1, r2, r3);
                   double N_temp = 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                   if (r1 != r2 && r1 != r3 && r2 != r3) {
                       V = npcf::gaussQuadCrossSectionDDR(DD, r2, r3, r1);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(DD, r3, r1, r2);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(DD, r1, r3, r2);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(DD, r2, r1, r3);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(DD, r3, r2, r1);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                   } else if ((r1 == r2 && r1 != r3) || (r1 == r3 && r1 != r2) || (r2 == r3 && r2 != r1)) {
                       V = npcf::gaussQuadCrossSectionDDR(DD, r2, r3, r1);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                       V = npcf::gaussQuadCrossSectionDDR(DD, r3, r1, r2);
                       N_temp += 4.0*M_PI*npcf::nbar_ran*V*npcf::N_parts;
                   }
                   npcf::DDR[index] = int(floor(N_temp + 0.5));
                }
            }
        }
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

int npcf::calculateCorrelations(float3 *galaxies) {
    float l = (float)pow(npcf::V_box, 1.0/3.0);
    float3 L = {l, l, l};
    float R = (float)npcf::r_max;
    float r = (float)npcf::r_min;
    std::vector<int3> shifts = getShifts();
    
    cudaMemcpyToSymbol(d_shift, shifts.data(), shifts.size()*sizeof(int3));
    cudaMemcpyToSymbol(d_L, &L, sizeof(float3));
    cudaMemcpyToSymbol(d_R, &R, sizeof(float));
    cudaMemcpyToSymbol(d_r, &r, sizeof(float));
    cudaMemcpyToSymbol(d_Nparts, &npcf::N_parts, sizeof(int));
    cudaMemcpyToSymbol(d_Nshells, &npcf::N_shells, sizeof(int));
    
    // Find the cell size such that the box size is an integer multiple of the cell size
    int l_c = floor(L.x/R);
    double L_c = L.x/l_c;
    int3 N = {l_c, l_c, l_c};
