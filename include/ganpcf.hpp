#ifndef _GANPCF_HPP_
#define _GANPCF_HPP_

#include <vector>
#include <vector_types.h>

class npcf{
    // Private data members
    double nbar_ran, Delta_r, r_max, r_min, V_box;
    int N_parts, N_rans, N_shells, timesRan;
    std::vector<double> rs, w, x;
    std::vector<double> twoPoint, threePoint;
    std::vector<float3> triangles;
    std::vector<int> DD, DR, DDD, DDR, DRR, RRR;
    
    
    // Private member functions
    void initVectors();
    
    void rezeroVectors();
    
    void swapIfGreater(double &a, double &b);
    
    double sphereOverlapVolume(double d, double R, double r);
    
    double crossSectionVolume(double r1, double r2, double r3);
    
    int getPermutations(double r1, double r2, double r3);
    
    double sphericalShellVolume(double r);
    
    double nbarData(std::vector<int> &DD, double r, double r1);
    
    double gaussQuadCrossSection(double r1, double r2, double r3);
    
    double gaussQuadCrossSectionDDR(std::vector<int> &DD, double r1, double r2, double r3);
    
    void getRRR();
    
    void getDRR();
    
    void getDDR();
    
    void getDR();
    
    // Public class interface
    public:
        npcf(int timesRandoms, int numShells, double volBox, double rMin, double rMax);
        
        int getShells(double *shells[]);
        
        int getNumTriangles();
        
        int getTriangles(float3 *tris[]);
        
        int setNumParticles(int numParticles);
        
        int get2ptSize();
        
        int get3ptSize();
        
        int calculateCorrelations(float3 *galaxies[]);
        
        int get2pt(double *twoPt[]);
        
        int get3pt(double *threePt[]);
};

#endif
