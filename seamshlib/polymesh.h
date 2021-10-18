#ifndef POLYMESH_H
#define POLYMESH_H
typedef struct PolyMesh PolyMesh;
PolyMesh *polymesh_new(double xmin[2], double xmax[2]);
void polymesh_add_points(PolyMesh *pm, int n, double *x, int *tags);
void polymesh_delete(PolyMesh *pm);
int polymesh_n_faces(const PolyMesh *pm);
void polymesh_faces(const PolyMesh *pm, int *faces);
#endif
