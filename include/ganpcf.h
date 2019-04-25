#ifndef _GANPCF_H_
#define _GANPCF_H_

#ifdef __cplusplus
extern "C" {
    struct float3;
    class npcf;
    typedef npcf NPCF;
#else
    typedef struct NPCF NPCF;
#endif
    
    NPCF* create_npcf(int timesRandoms, int numShells, double volBox, double rMin, double rMax);
    
    void delete_npcf(NPCF* obj);
    
    int get_shells(NPCF *obj, double *shells[]);
    
    int get_num_triangles(NPCF *obj);
    
    int get_triangles(NPCF *obj, float3 *tris[]);
    
    int set_num_particles(NPCF *obj, int numParticles);
    
    int get_2pt_size(NPCF *obj);
    
    int get_3pt_size(NPCF *obj);
    
    int calculate_correlations(NPCF *obj, float3 *galaxies[]);
    
    int calculate_2pt(NPCF *obj, float3 *galaxies[]);
    
    int get_2pt(NPCF *obj, double *twoPt[]);
    
    int get_3pt(NPCF *obj, double *threePt[]);
    
#ifdef __cplusplus
}
#endif

#endif
