// Gmsh - Copyright (C) 1997-2021 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include "vector.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "robustPredicates.h"
#include "polymesh.h"

typedef struct HalfEdge HalfEdge;
typedef struct Face Face;
typedef struct Vertex Vertex;

struct Vertex{
  double p[2];
  int data;
  HalfEdge *he; // one incident half edge
};

struct HalfEdge{
  Vertex *v; // origin
  Face *f; // incident face
  HalfEdge *prev; // previous half edge on the face
  HalfEdge *next; // next half edge on the face
  HalfEdge *opposite; // opposite half edge (twin)
  int data;
};

struct Face {
  HalfEdge *he;
  int data;
};

struct PolyMesh {
  Vertex **vertices;
  HalfEdge **hedges;
  Face **faces;
};

Vertex *vertex_new(double x, double y, int d)
{
  Vertex *v = (Vertex*)malloc(sizeof(Vertex));
  v->p[0] = x;
  v->p[1] = y;
  v->data = d;
  v->he = NULL;
  return v;
}


HalfEdge *half_edge_new(Vertex *v) {
  HalfEdge *he = (HalfEdge*) malloc(sizeof(HalfEdge));
  he->v = v;
  he->f = NULL;
  he->prev = NULL;
  he->next = NULL;
  he->opposite = NULL;
  he->data = -1;
  return he;
}

void half_edge_dir(const HalfEdge *he, double dir[2]) {
  dir[0] = he->next->v->p[0] - he->v->p[1];
  dir[1] = he->next->v->p[1] - he->v->p[1];
  double l = hypot(dir[0], dir[1]);
  dir[0] /= l;
  dir[1] /= l;
}

Face *face_new(HalfEdge *e) {
  Face *f = (Face*)malloc(sizeof(Face));
  f->he = e;
  f->data = -1;
  return f;
}


PolyMesh *_poly_mesh_new() {
  PolyMesh *pm = (PolyMesh*)malloc(sizeof(PolyMesh));
  pm->vertices = NULL;
  pm->hedges = NULL;
  pm->faces = NULL;
  return pm;
}

void poly_mesh_reset(PolyMesh *pm)
{
  for (size_t i = 0; i < vector_size(pm->vertices); ++i) {
    free((void*)pm->vertices[i]);
  }
  vector_free(pm->vertices);
  for (size_t i = 0; i < vector_size(pm->hedges); ++i) {
    free((void*)pm->hedges[i]);
  }
  vector_free(pm->hedges);
  for (size_t i = 0; i < vector_size(pm->faces); ++i) {
    free((void*)pm->faces[i]);
  }
  vector_free(pm->faces);
}

void poly_mesh_delete(PolyMesh *pm) {
  poly_mesh_reset(pm);
  free(pm);
}

// compute the degree of a given vertex v
static inline int degree(const Vertex *v)
{
  HalfEdge *he = v->he;
  size_t count = 0;
  do {
    he = he->opposite;
    if(he == NULL) return -1;
    he = he->next;
    count++;
  } while(he != v->he);
  return count;
}

static inline int num_sides(const HalfEdge *he)
{
  size_t count = 0;
  const HalfEdge *start = he;
  do {
    count++;
    he = he->next;
  } while(he != start);
  return count;
}

static inline HalfEdge *getEdge(Vertex *v0, Vertex *v1)
{
  HalfEdge *he = v0->he;
  do {
    if(he->next->v == v1) return he;
    he = he->opposite;
    if(he == NULL) return NULL;
    he = he->next;
  } while(he != v0->he);
  return NULL;
}

static inline void createFace(Face *f, Vertex *v0, Vertex *v1, Vertex *v2,
    HalfEdge *he0, HalfEdge *he1, HalfEdge *he2)
{
  he0->v = v0;
  he1->v = v1;
  he2->v = v2;
  v0->he = he0;
  v1->he = he1;
  v2->he = he2;

  he0->next = he1;
  he1->prev = he0;
  he1->next = he2;
  he2->prev = he1;
  he2->next = he0;
  he0->prev = he2;
  he0->f = he1->f = he2->f = f;
  f->he = he0;
}

// swap without asking questions
//
//         he1
// v2 ------->------ v3
//    | \          |
//    |   \ he0    | he2
// heo2|     \      |
//    |  heo0 \    |
//    |         \  |
// v1 -----<------- v0
//          heo1
//
//          he1
//    --------------
//    |         /  |
//    |   he0 /    | he2
// heo2|    /       |
//    |  /heo0     |
//    |/           |
//    --------------
//          heo1
//

static inline int swap_edge(HalfEdge *he0)
{
  HalfEdge *heo0 = he0->opposite;
  if(heo0 == NULL) return -1;

  HalfEdge *he1 = he0->next;
  HalfEdge *he2 = he1->next;
  HalfEdge *heo1 = heo0->next;
  HalfEdge *heo2 = heo1->next;

  Vertex *v0 = heo1->v;
  Vertex *v1 = heo2->v;
  Vertex *v2 = heo0->v;
  Vertex *v3 = he2->v;

  createFace(he0->f, v0, v1, v3, heo1, heo0, he2);
  createFace(heo2->f, v1, v2, v3, heo2, he1, he0);
  return 0;
}

static inline int merge_faces(HalfEdge *he)
{
  HalfEdge *heo = he->opposite;

  if(heo == NULL) return -1;

  Face *to_delete = heo->f;

  do {
    heo->f = he->f;
    heo = heo->next;
  } while(heo != he->opposite);

  he->next->prev = heo->prev;
  heo->prev->next = he->next;
  he->prev->next = heo->next;
  heo->next->prev = he->prev;

  he->f->he = he->next;
  he->v->he = heo->next;
  heo->v->he = he->next;

  // remove afterwards...
  he->v = NULL;
  heo->v = NULL;
  to_delete->he = NULL;
  return 0;
}

static void cleanv(PolyMesh *pm)
{
  Vertex **uv = NULL;
  for(size_t i = 0; i < vector_size(pm->vertices); ++i) {
    Vertex *v = pm->vertices[i];
    if(v->he)
      *vector_push(&uv) = v;
    else
      free(v);
  }
  vector_free(pm->vertices);
  pm->vertices = uv;
}

static void cleanh(PolyMesh *pm)
{
  HalfEdge **uh = NULL;
  for(size_t i = 0; i < vector_size(pm->hedges); ++i) {
    HalfEdge *h = pm->hedges[i];
    if(h->f)
      *vector_push(&uh) = h;
    else
      free(h);
  }
  vector_free(pm->hedges);
  pm->hedges = uh;
}

static void cleanf(PolyMesh *pm)
{
  Face **uf = NULL;
  for(size_t i = 0; i < vector_size(pm->faces); ++i) {
    Face *f = pm->faces[i];
    if(f->he)
      *vector_push(&uf) = f;
    else
      free(f);
  }
  vector_free(pm->faces);
  pm->faces = uf;
}

static void clean(PolyMesh *pm)
{
  cleanv(pm);
  cleanh(pm);
  cleanf(pm);
}

static inline int split_edge(PolyMesh *pm, HalfEdge *he0m, const double position[2], int data)
{
  HalfEdge *he1m = he0m->opposite;
  if(he1m == NULL) return -1;

  Vertex *mid = vertex_new(position[0], position[1], data);
  *vector_push(&pm->vertices) = mid;

  HalfEdge *he12 = he0m->next;
  HalfEdge *he20 = he0m->next->next;
  HalfEdge *he03 = he0m->opposite->next;
  HalfEdge *he31 = he0m->opposite->next->next;

  // if(he03->v != he0m->v) Msg::Error("error 1");
  // if(he1m->v != he12->v) Msg::Error("error 2");

  Vertex *v0 = he03->v;
  Vertex *v1 = he12->v;
  Vertex *v2 = he20->v;
  Vertex *v3 = he31->v;

  HalfEdge *hem0 = half_edge_new(mid);
  HalfEdge *hem1 = half_edge_new(mid);
  HalfEdge *hem2 = half_edge_new(mid);
  HalfEdge *hem3 = half_edge_new(mid);

  HalfEdge *he2m = half_edge_new(v2);
  HalfEdge *he3m = half_edge_new(v3);

  he0m->opposite = hem0;
  hem0->opposite = he0m;
  he1m->opposite = hem1;
  hem1->opposite = he1m;
  he2m->opposite = hem2;
  hem2->opposite = he2m;
  he3m->opposite = hem3;
  hem3->opposite = he3m;

  *vector_push(&pm->hedges) = hem0;
  *vector_push(&pm->hedges) = hem1;
  *vector_push(&pm->hedges) = hem2;
  *vector_push(&pm->hedges) = hem3;
  *vector_push(&pm->hedges) = he2m;
  *vector_push(&pm->hedges) = he3m;

  Face *f0m2 = he0m->f;
  Face *f1m3 = he1m->f;
  Face *f2m1 = face_new(he2m);
  Face *f3m0 = face_new(he3m);
  *vector_push(&pm->faces) = f2m1;
  *vector_push(&pm->faces) = f3m0;

  createFace(f0m2, v0, mid, v2, he0m, hem2, he20);
  createFace(f1m3, v1, mid, v3, he1m, hem3, he31);
  createFace(f2m1, v2, mid, v1, he2m, hem1, he12);
  createFace(f3m0, v3, mid, v0, he3m, hem0, he03);
  return 0;
}

//
// v0   he0
// ------------------>------ v1
// |                      /
// |                   /
// |      v         /
// |he2          /
// |          /  he1
// |       /
// |    /
// |/
// v2

static inline void initialize_rectangle(PolyMesh *pm, double xmin, double xmax, double ymin,
    double ymax)
{
  poly_mesh_reset(pm);
  Vertex *v_mm = vertex_new(xmin, ymin, -1);
  *vector_push(&pm->vertices) = v_mm;
  Vertex *v_mM = vertex_new(xmin, ymax, -1);
  *vector_push(&pm->vertices) = v_mM;
  Vertex *v_MM = vertex_new(xmax, ymax, -1);
  *vector_push(&pm->vertices) = v_MM;
  Vertex *v_Mm = vertex_new(xmax, ymin, -1);
  *vector_push(&pm->vertices) = v_Mm;
  HalfEdge *mm_MM = half_edge_new(v_mm);
  HalfEdge *MM_Mm = half_edge_new(v_MM);
  HalfEdge *Mm_mm = half_edge_new(v_Mm);
  *vector_push(&pm->hedges) = mm_MM;
  *vector_push(&pm->hedges) = MM_Mm;
  *vector_push(&pm->hedges) = Mm_mm;
  Face *f0 = face_new(mm_MM);
  *vector_push(&pm->faces) = f0;
  createFace(f0, v_mm, v_MM, v_Mm, mm_MM, MM_Mm, Mm_mm);

  HalfEdge *MM_mm = half_edge_new(v_MM);
  HalfEdge *mm_mM = half_edge_new(v_mm);
  HalfEdge *mM_MM = half_edge_new(v_mM);
  *vector_push(&pm->hedges) = MM_mm;
  *vector_push(&pm->hedges) = mm_mM;
  *vector_push(&pm->hedges) = mM_MM;
  Face *f1 = face_new(MM_mm);
  *vector_push(&pm->faces) = f1;
  createFace(f1, v_MM, v_mm, v_mM, MM_mm, mm_mM, mM_MM);

  MM_mm->opposite = mm_MM;
  mm_MM->opposite = MM_mm;
}

static inline int split_triangle(PolyMesh *pm, int index, double x, double y, Face *f,
    int (*doSwap)(HalfEdge *, void *),
    void *data,
    HalfEdge ***_t)
{
  Vertex *v = vertex_new(x, y, -1); // one more vertex
  *vector_push(&pm->vertices) = v;

  HalfEdge *he0 = f->he;
  HalfEdge *he1 = he0->next;
  HalfEdge *he2 = he1->next;

  Vertex *v0 = he0->v;
  Vertex *v1 = he1->v;
  Vertex *v2 = he2->v;
  HalfEdge *hev0 = half_edge_new(v);
  HalfEdge *hev1 = half_edge_new(v);
  HalfEdge *hev2 = half_edge_new(v);

  HalfEdge *he0v = half_edge_new(v0);
  HalfEdge *he1v = half_edge_new(v1);
  HalfEdge *he2v = half_edge_new(v2);

  *vector_push(&pm->hedges) = hev0;
  *vector_push(&pm->hedges) = hev1;
  *vector_push(&pm->hedges) = hev2;
  *vector_push(&pm->hedges) = he0v;
  *vector_push(&pm->hedges) = he1v;
  *vector_push(&pm->hedges) = he2v;

  hev0->opposite = he0v;
  he0v->opposite = hev0;
  hev1->opposite = he1v;
  he1v->opposite = hev1;
  hev2->opposite = he2v;
  he2v->opposite = hev2;

  Face *f0 = f;
  f->he = hev0;
  Face *f1 = face_new(hev1);
  Face *f2 = face_new(hev2);
  f1->data = f2->data = f0->data;

  *vector_push(&pm->faces) = f1;
  *vector_push(&pm->faces) = f2;

  createFace(f0, v0, v1, v, he0, he1v, hev0);
  createFace(f1, v1, v2, v, he1, he2v, hev1);
  createFace(f2, v2, v0, v, he2, he0v, hev2);

  if(doSwap) {
    HalfEdge **_stack = NULL;
    *vector_push(&_stack) = he0;
    *vector_push(&_stack) = he1;
    *vector_push(&_stack) = he2;
    HalfEdge **touched = NULL;
    while(vector_size(_stack) != 0) {
      HalfEdge *he = _stack[vector_size(_stack)-1];
      vector_pop(_stack);
      *vector_push(&touched) = he;
      //	printf("do we swap %g %g --> %g %g ?\n",
      //		       he->v->position.x(),he->v->position.y(),
      //	he->next->v->position.x(),he->next->v->position.y());
      if(doSwap(he, data) == 1) {
        //	  printf("YES\n");
        swap_edge(he);

        HalfEdge *H[2] = {he, he->opposite};

        for(int k = 0; k < 2; k++) {
          if(H[k] == NULL) continue;
          HalfEdge *heb = H[k]->next;
          HalfEdge *hebo = heb->opposite;
          HalfEdge *hec = heb->next;
          HalfEdge *heco = hec->opposite;

          int found_b_bo = 0;
          int found_c_co = 0;
          for (size_t i = 0; i < vector_size(touched); ++i) {
            if (touched[i] == heb || touched[i] == hebo)
              found_b_bo = 1;
            if (touched[i] == hec || touched[i] == heco)
              found_c_co = 1;
          }

          if(found_b_bo == 0) {
            *vector_push(&_stack) = heb;
          }
          if (found_c_co == 0) {
            *vector_push(&_stack) = hec;
          }
        }
      }
    }
    if(_t) *_t = touched;
    else vector_free(touched);
    vector_free(_stack);
  }
  return 0;
}

inline Vertex *ptrmin(Vertex *p1, Vertex *p2) {
  return p1 < p2 ? p1 : p2;
}

inline Vertex *ptrmax(Vertex *p1, Vertex *p2) {
  return p1 > p2 ? p1 : p2;
}
/*
struct HalfEdgePtrLessThan {
  bool operator()(HalfEdge *l1, HalfEdge *l2) const
  {
    Vertex *l10 = ptrmin(l1->v, l1->next->v);
    Vertex *l11 = ptrmax(l1->v, l1->next->v);
    Vertex *l20 = ptrmin(l2->v, l2->next->v);
    Vertex *l21 = ptrmax(l2->v, l2->next->v);
    if(l10 < l20) return true;
    if(l10 > l20) return false;
    if(l11 > l21) return true;
    return false;
  }
};

struct HalfEdgePtrEqual {
  bool operator()(HalfEdge *l1, HalfEdge *l2) const
  {
    Vertex *l10 = ptrmin(l1->v, l1->next->v);
    Vertex *l11 = ptrmax(l1->v, l1->next->v);
    Vertex *l20 = ptrmin(l2->v, l2->next->v);
    Vertex *l21 = ptrmax(l2->v, l2->next->v);
    if(l10 == l20 && l11 == l21) return true;
    return false;
  }
};*/


void swap(double *a, double *b)
{
  double temp = *a;
  *a = *b;
  *b = temp;
}

size_t HilbertCoordinates(double x, double y, double x0, double y0, double xRed,
                          double yRed, double xBlue, double yBlue)
{
  size_t BIG = 1073741824;
  size_t RESULT = 0;
  for(int i = 0; i < 16; i++) {
    double coordRed = (x - x0) * xRed + (y - y0) * yRed;
    double coordBlue = (x - x0) * xBlue + (y - y0) * yBlue;
    xRed /= 2;
    yRed /= 2;
    xBlue /= 2;
    yBlue /= 2;
    if(coordRed <= 0 && coordBlue <= 0) { // quadrant 0
      x0 -= (xBlue + xRed);
      y0 -= (yBlue + yRed);
      swap(&xRed, &xBlue);
      swap(&yRed, &yBlue);
    }
    else if(coordRed <= 0 && coordBlue >= 0) { // quadrant 1
      RESULT += BIG;
      x0 += (xBlue - xRed);
      y0 += (yBlue - yRed);
    }
    else if(coordRed >= 0 && coordBlue >= 0) { // quadrant 2
      RESULT += 2 * BIG;
      x0 += (xBlue + xRed);
      y0 += (yBlue + yRed);
    }
    else if(coordRed >= 0 && coordBlue <= 0) { // quadrant 3
      x0 += (-xBlue + xRed);
      y0 += (-yBlue + yRed);
      swap(&xRed, &xBlue);
      swap(&yRed, &yBlue);
      xBlue = -xBlue;
      yBlue = -yBlue;
      xRed = -xRed;
      yRed = -yRed;
      RESULT += 3 * BIG;
    }
    else
      printf("Hilbert failed %g %g", coordRed, coordBlue);
    BIG /= 4;
  }
  return RESULT;
}

static Face *Walk(Face *f, double x, double y)
{
  double POS[2] = {x, y};
  HalfEdge *he = f->he;

  while(1) {
    Vertex *v0 = he->v;
    Vertex *v1 = he->next->v;
    Vertex *v2 = he->next->next->v;

    double s0 = -orient2d(v0->p, v1->p, POS);
    double s1 = -orient2d(v1->p, v2->p, POS);
    double s2 = -orient2d(v2->p, v0->p, POS);

    if(s0 >= 0 && s1 >= 0 && s2 >= 0) {
      return he->f;
    }
    else if(s0 <= 0 && s1 >= 0 && s2 >= 0)
      he = he->opposite;
    else if(s1 <= 0 && s0 >= 0 && s2 >= 0)
      he = he->next->opposite;
    else if(s2 <= 0 && s0 >= 0 && s1 >= 0)
      he = he->next->next->opposite;
    else if(s0 <= 0 && s1 <= 0)
      he = s0 > s1 ? he->opposite : he->next->opposite;
    else if(s0 <= 0 && s2 <= 0)
      he = s0 > s2 ? he->opposite : he->next->next->opposite;
    else if(s1 <= 0 && s2 <= 0)
      he = s1 > s2 ? he->next->opposite : he->next->next->opposite;
    else {
      printf("Could not find half-edge in walk for point %g %g on "
                 "face %g %g / %g %g / %g %g "
                 "(orientation tests %g %g %g)", x, y,
                 v0->p[0], v0->p[1],
                 v1->p[0], v1->p[1],
                 v2->p[0], v2->p[1],
                 s0, s1, s2);
    }
    if(he == NULL) break;
  }
  // should only come here wether the triangulated domain is not convex
  return NULL;
}

static int delaunayEdgeCriterionPlaneIsotropic(HalfEdge *he, void *ignore)
{
  if(he->opposite == NULL) return -1;
  Vertex *v0 = he->v;
  Vertex *v1 = he->next->v;
  Vertex *v2 = he->next->next->v;
  Vertex *v = he->opposite->next->next->v;

  // FIXME : should be oriented anyway !
  double result = -incircle(v0->p, v1->p, v2->p, v->p);

  return (result > 0) ? 1 : 0;
}

int polymesh_n_faces(const PolyMesh *pm) {
  int n = 0;
  for (size_t i = 0; i < vector_size(pm->faces); ++i) {
    HalfEdge *he = pm->faces[i]->he;
    int tri[3];
    for (int j = 0; j < 3; ++j) {
      tri[j] = he->v->data;
      he = he->next;
    }
    n += (tri[0] >= 0 && tri[1] >= 0 && tri[2] >= 0);
  }
  return n;
}

void polymesh_faces(const PolyMesh *pm, int *faces) {
  int n = 0;
  for (size_t i = 0; i < vector_size(pm->faces); ++i) {
    HalfEdge *he = pm->faces[i]->he;
    int tri[3];
    for (int j = 0; j < 3; ++j) {
      tri[j] = he->v->data;
      he = he->next;
    }
    if (tri[0] >= 0 && tri[1] >= 0 && tri[2] >= 0) {
      faces[n*3+0] = tri[0];
      faces[n*3+1] = tri[1];
      faces[n*3+2] = tri[2];
      n++;
    }
  }
}

void polymesh_delete(PolyMesh *pm) {
  poly_mesh_reset(pm);
  free(pm);
}

static void get_bounding_box(int n, double *x, double bbmin[2], double bbmax[2]) {
  bbmin[0] = bbmin[1] = DBL_MAX;
  bbmax[0] = bbmax[1] = -DBL_MAX;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < 2; ++j) {
      bbmin[j] = fmin(x[i*2+j], bbmin[j]);
      bbmax[j] = fmax(x[i*2+j], bbmax[j]);
    }
  }
  for (int j = 0; j < 2; ++j) {
    double L = bbmax[j]-bbmin[j];
    bbmin[j] -= L*0.1;
    bbmax[j] += L*0.1;
  }
}

PolyMesh *polymesh_new(double xmin[2], double xmax[2]) {
  PolyMesh *pm = _poly_mesh_new();
  double d[2] = {xmax[0]-xmin[0],xmax[1]-xmin[1]};
  initialize_rectangle(pm, xmin[0]-d[0]*0.1, xmax[0]+d[0]*0.1,
      xmin[1]-d[1]*0.1, xmax[1]+d[1]*0.1);
  return pm;
}


void polymesh_add_points(PolyMesh *pm, int n, double *x, int *tags)
{
  size_t *HC = (size_t*)malloc(sizeof(size_t)*n);
  size_t *IND = NULL;
  vector_push_n(&IND, n);
  Face *f = pm->faces[0];
  double bbmin[2], bbmax[2];
  get_bounding_box(n, x, bbmin, bbmax);
  double bbcenter[2] = {(bbmin[0]+bbmax[0])/2, (bbmin[1]+bbmax[1])/2};
  for(int i = 0; i < n; i++) {
    HC[i] = HilbertCoordinates(x[i*2], x[i*2+1], bbcenter[0], bbcenter[1],
                               bbmax[0] - bbcenter[0], 0, 0,
                               bbmax[1] - bbcenter[1]);
    IND[i] = i;
  }
  #define cmp(i,j) (HC[IND[i]]<HC[IND[j]])
  vector_sort_r(IND, cmp);
  for(int i = 0; i < n; i++) {
    size_t I = IND[i];
    f = Walk(f, x[I*2], x[I*2+1]);
    split_triangle(pm, i, x[I*2], x[I*2+1], f, delaunayEdgeCriterionPlaneIsotropic, NULL, NULL);
    pm->vertices[vector_size(pm->vertices) - 1]->data = tags[I];
  }
  free(HC);
  vector_free(IND);
}
