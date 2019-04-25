#include "../include/ganpcf.hpp"
#include "../include/ganpcf.h"

NPCF *create_npcf(int timesRandoms, int numShells, double volBox, double rMin, double rMax) {
    return new npcf(timesRandoms, numShells, volBox, rMin, rMax);
}

void delete_npcf(NPCF *obj) {
    delete obj;
}

int get_shells(NPCF *obj, double *shells[]) {
    return obj->getShells(shells);
}

int get_num_triangles(NPCF *obj) {
    return obj->getNumTriangles();
}

int get_triangles(NPCF *obj, float3 *tris[]) {
    return obj->getTriangles(tris);
}

int set_num_particles(NPCF *obj, int numParticles) {
    return obj->setNumParticles(numParticles);
}

int get_2pt_size(NPCF *obj) {
    return obj->get2ptSize();
}

int get_3pt_size(NPCF *obj) {
    return obj->get3ptSize();
}

int calculate_correlations(NPCF *obj, float3 *galaxies[]) {
    return obj->calculateCorrelations(galaxies);
}

int calculate_2pt(NPCF *obj, float3 *galaxies[]) {
    return obj->calculate2pt(galaxies);
}

int get_2pt(NPCF *obj, double *twoPt[]) {
    return obj->get2pt(twoPt);
}

int get_3pt(NPCF *obj, double *threePt[]) {
    return obj->get3pt(threePt);
}
