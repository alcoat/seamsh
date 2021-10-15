// Gmsh - Copyright (C) 1997-2021 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <vector>
#include <algorithm>
#include <stack>
#include <stdio.h>
#include <cmath>
#include <limits>
#include "robustPredicates.h"

class SVector3 {
  double _x[3];
  public:
  SVector3(double x, double y, double z){
    _x[0] = x;
    _x[1] = y;
    _x[2] = z;
  }
  SVector3(const SVector3 &o) {
    _x[0] = o._x[0];
    _x[1] = o._x[1];
    _x[2] = o._x[2];
  }
  SVector3 operator-(const SVector3 &o) const {
    return SVector3(_x[0]-o._x[0],_x[1]-o._x[1],_x[2]-o._x[2]);
  }
  SVector3 operator+(const SVector3 &o) const {
    return SVector3(_x[0]+o._x[0],_x[1]+o._x[1],_x[2]+o._x[2]);
  }
  void operator +=(const SVector3 &o) {
    for (int i=0; i<3; ++i) _x[i] += o._x[i];
  }
  SVector3 (double v) {
    _x[0] = _x[1] = _x[2] = v;
  }
  SVector3 operator/(double d) const {
    return SVector3(x()/d,y()/d,z()/d);
  }
  operator double*(){return _x;}
  void normalize() {
    double n = sqrt(_x[0]*_x[0]+_x[1]*_x[1]+_x[2]*_x[2]);
    _x[0] /= n;
    _x[1] /= n;
    _x[2] /= n;
  }
  double x()const {return _x[0];}
  double y()const {return _x[1];}
  double z()const {return _x[2];}
};

SVector3 crossprod(const SVector3 &a, const SVector3 &b) {
  return SVector3( a.y()*b.z() - a.z()*b.y(),
                   a.z()*b.x() - a.x()*b.z(),
                   a.x()*b.y() - a.y()*b.x());
}

class PolyMesh {
public:
  class HalfEdge;
  class Face;
  class Vertex;

  class Vertex {
  public:
    Vertex(double x, double y, double z, int _d = -1)
      : position(x, y, z), he(NULL), data(_d)
    {
    }
    SVector3 position;
    PolyMesh::HalfEdge *he; // one incident half edge
    int data;
  };

  class HalfEdge {
  public:
    HalfEdge(Vertex *vv)
      : v(vv), f(NULL), prev(NULL), next(NULL), opposite(NULL), data(-1)
    {
    }
    Vertex *v; // origin
    Face *f; // incident face
    HalfEdge *prev; // previous half edge on the face
    HalfEdge *next; // next half edge on the face
    HalfEdge *opposite; // opposite half edge (twin)
    int data;
    SVector3 d() const
    {
      SVector3 t = next->v->position - v->position;
      t.normalize();
      return t;
    }
  };

  class Face {
  public:
    Face(HalfEdge *e) : he(e), data(-1) {}
    HalfEdge *he; // one half edge of the face
    int data;
  };

  std::vector<Vertex *> vertices;
  std::vector<HalfEdge *> hedges;
  std::vector<Face *> faces;
  std::vector<SVector3> high_order_nodes;

  void reset()
  {
    for(auto it : vertices) delete it;
    for(auto it : hedges) delete it;
    for(auto it : faces) delete it;
  }

  ~PolyMesh() { reset(); }

  void print4debug(const int debugTag)
  {
    char name[256];
    sprintf(name, "polyMesh%d.pos", debugTag);
    FILE *f = fopen(name, "w");
    fprintf(f, "View \" %s \"{\n", name);
    for(auto it : faces) {
      HalfEdge *he0 = it->he;
      HalfEdge *he1 = it->he->next;
      HalfEdge *he2 = it->he->next->next;
      fprintf(f, "ST(%g,%g,0,%g,%g,0,%g,%g,0){%d,%d,%d};\n",
              he0->v->position.x(), he0->v->position.y(), he1->v->position.x(),
              he1->v->position.y(), he2->v->position.x(), he2->v->position.y(),
              it->data, it->data, it->data);
    }
    for(auto it : hedges) {
      HalfEdge *he = it;
      if(he->data >= 0) {
        fprintf(f, "SL(%g,%g,0,%g,%g,0){%d,%d};\n", he->v->position.x(),
                he->v->position.y(), he->opposite->v->position.x(),
                he->opposite->v->position.y(), he->data, he->data);
      }
    }

    fprintf(f, "};\n");
    fclose(f);
  }

  // compute the degree of a given vertex v
  inline int degree(const Vertex *v) const
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

  inline int num_sides(const HalfEdge *he) const
  {
    size_t count = 0;
    const HalfEdge *start = he;
    do {
      count++;
      he = he->next;
    } while(he != start);
    return count;
  }

  // compute the normal of an internal vertex v
  inline SVector3 normal(const Vertex *v) const
  {
    SVector3 n(0, 0, 0);
    HalfEdge *he = v->he;
    do {
      SVector3 n1 = he->d();
      he = he->opposite;
      if(he == NULL) return -1;
      he = he->next;
      n += crossprod(n1, he->d());
    } while(he != v->he);
    n.normalize();
    return n;
  }

  inline HalfEdge *getEdge(Vertex *v0, Vertex *v1)
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

  inline void createFace(Face *f, Vertex *v0, Vertex *v1, Vertex *v2,
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

  inline int swap_edge(HalfEdge *he0)
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

  inline int merge_faces(HalfEdge *he)
  {
    PolyMesh::HalfEdge *heo = he->opposite;

    if(heo == nullptr) return -1;

    PolyMesh::Face *to_delete = heo->f;

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
    he->v = nullptr;
    heo->v = nullptr;
    to_delete->he = nullptr;
    return 0;
  }

  void cleanv()
  {
    std::vector<Vertex *> uv;
    for(auto v : vertices) {
      if(v->he)
        uv.push_back(v);
      else
        delete v;
    }
    vertices = uv;
  }

  void cleanh()
  {
    std::vector<HalfEdge *> uh;
    for(auto h : hedges) {
      if(h->f)
        uh.push_back(h);
      else
        delete h;
    }
    hedges = uh;
  }

  void cleanf()
  {
    std::vector<Face *> uf;
    for(auto f : faces) {
      if(f->he)
        uf.push_back(f);
      else
        delete f;
    }
    faces = uf;
  }

  void clean()
  {
    cleanv();
    cleanh();
    cleanf();
  }

  inline int split_edge(HalfEdge *he0m, const SVector3 &position, int data)
  {
    HalfEdge *he1m = he0m->opposite;
    if(he1m == nullptr) return -1;

    Vertex *mid = new Vertex(position.x(), position.y(), position.z(), data);
    vertices.push_back(mid);

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

    HalfEdge *hem0 = new HalfEdge(mid);
    HalfEdge *hem1 = new HalfEdge(mid);
    HalfEdge *hem2 = new HalfEdge(mid);
    HalfEdge *hem3 = new HalfEdge(mid);

    HalfEdge *he2m = new HalfEdge(v2);
    HalfEdge *he3m = new HalfEdge(v3);

    he0m->opposite = hem0;
    hem0->opposite = he0m;
    he1m->opposite = hem1;
    hem1->opposite = he1m;
    he2m->opposite = hem2;
    hem2->opposite = he2m;
    he3m->opposite = hem3;
    hem3->opposite = he3m;

    hedges.push_back(hem0);
    hedges.push_back(hem1);
    hedges.push_back(hem2);
    hedges.push_back(hem3);
    hedges.push_back(he2m);
    hedges.push_back(he3m);

    Face *f0m2 = he0m->f;
    Face *f1m3 = he1m->f;
    Face *f2m1 = new Face(he2m);
    Face *f3m0 = new Face(he3m);
    faces.push_back(f2m1);
    faces.push_back(f3m0);

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

  inline void initialize_rectangle(double xmin, double xmax, double ymin,
                                   double ymax)
  {
    reset();
    Vertex *v_mm = new PolyMesh::Vertex(xmin, ymin, 0);
    vertices.push_back(v_mm);
    Vertex *v_mM = new PolyMesh::Vertex(xmin, ymax, 0);
    vertices.push_back(v_mM);
    Vertex *v_MM = new PolyMesh::Vertex(xmax, ymax, 0);
    vertices.push_back(v_MM);
    Vertex *v_Mm = new PolyMesh::Vertex(xmax, ymin, 0);
    vertices.push_back(v_Mm);
    HalfEdge *mm_MM = new HalfEdge(v_mm);
    HalfEdge *MM_Mm = new HalfEdge(v_MM);
    HalfEdge *Mm_mm = new HalfEdge(v_Mm);
    hedges.push_back(mm_MM);
    hedges.push_back(MM_Mm);
    hedges.push_back(Mm_mm);
    Face *f0 = new Face(mm_MM);
    faces.push_back(f0);
    createFace(f0, v_mm, v_MM, v_Mm, mm_MM, MM_Mm, Mm_mm);

    HalfEdge *MM_mm = new HalfEdge(v_MM);
    HalfEdge *mm_mM = new HalfEdge(v_mm);
    HalfEdge *mM_MM = new HalfEdge(v_mM);
    hedges.push_back(MM_mm);
    hedges.push_back(mm_mM);
    hedges.push_back(mM_MM);
    Face *f1 = new Face(MM_mm);
    faces.push_back(f1);
    createFace(f1, v_MM, v_mm, v_mM, MM_mm, mm_mM, mM_MM);

    MM_mm->opposite = mm_MM;
    mm_MM->opposite = MM_mm;
  }

  inline int split_triangle(int index, double x, double y, double z, Face *f,
                            int (*doSwap)(PolyMesh::HalfEdge *, void *) = NULL,
                            void *data = NULL,
                            std::vector<HalfEdge *> *_t = NULL)
  {
    Vertex *v = new PolyMesh::Vertex(x, y, z); // one more vertex
    v->data = -1;

    vertices.push_back(v);

    HalfEdge *he0 = f->he;
    HalfEdge *he1 = he0->next;
    HalfEdge *he2 = he1->next;

    Vertex *v0 = he0->v;
    Vertex *v1 = he1->v;
    Vertex *v2 = he2->v;
    HalfEdge *hev0 = new PolyMesh::HalfEdge(v);
    HalfEdge *hev1 = new PolyMesh::HalfEdge(v);
    HalfEdge *hev2 = new PolyMesh::HalfEdge(v);

    HalfEdge *he0v = new PolyMesh::HalfEdge(v0);
    HalfEdge *he1v = new PolyMesh::HalfEdge(v1);
    HalfEdge *he2v = new PolyMesh::HalfEdge(v2);

    hedges.push_back(hev0);
    hedges.push_back(hev1);
    hedges.push_back(hev2);
    hedges.push_back(he0v);
    hedges.push_back(he1v);
    hedges.push_back(he2v);

    hev0->opposite = he0v;
    he0v->opposite = hev0;
    hev1->opposite = he1v;
    he1v->opposite = hev1;
    hev2->opposite = he2v;
    he2v->opposite = hev2;

    Face *f0 = f;
    f->he = hev0;
    Face *f1 = new Face(hev1);
    Face *f2 = new Face(hev2);
    f1->data = f2->data = f0->data;

    faces.push_back(f1);
    faces.push_back(f2);

    createFace(f0, v0, v1, v, he0, he1v, hev0);
    createFace(f1, v1, v2, v, he1, he2v, hev1);
    createFace(f2, v2, v0, v, he2, he0v, hev2);

    if(doSwap) {
      std::stack<HalfEdge *> _stack;
      _stack.push(he0);
      _stack.push(he1);
      _stack.push(he2);
      std::vector<HalfEdge *> _touched;
      while(!_stack.empty()) {
        HalfEdge *he = _stack.top();
        _touched.push_back(he);
        _stack.pop();
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

            if(std::find(_touched.begin(), _touched.end(), heb) ==
                 _touched.end() &&
               std::find(_touched.begin(), _touched.end(), hebo) ==
                 _touched.end()) {
              _stack.push(heb);
            }

            HalfEdge *hec = heb->next;
            HalfEdge *heco = hec->opposite;

            if(std::find(_touched.begin(), _touched.end(), hec) ==
                 _touched.end() &&
               std::find(_touched.begin(), _touched.end(), heco) ==
                 _touched.end()) {
              _stack.push(hec);
            }
          }
        }
      }
      if(_t) *_t = _touched;
    }
    return 0;
  }
};

struct HalfEdgePtrLessThan {
  bool operator()(PolyMesh::HalfEdge *l1, PolyMesh::HalfEdge *l2) const
  {
    PolyMesh::Vertex *l10 = std::min(l1->v, l1->next->v);
    PolyMesh::Vertex *l11 = std::max(l1->v, l1->next->v);
    PolyMesh::Vertex *l20 = std::min(l2->v, l2->next->v);
    PolyMesh::Vertex *l21 = std::max(l2->v, l2->next->v);
    if(l10 < l20) return true;
    if(l10 > l20) return false;
    if(l11 > l21) return true;
    return false;
  }
};

struct HalfEdgePtrEqual {
  bool operator()(PolyMesh::HalfEdge *l1, PolyMesh::HalfEdge *l2) const
  {
    PolyMesh::Vertex *l10 = std::min(l1->v, l1->next->v);
    PolyMesh::Vertex *l11 = std::max(l1->v, l1->next->v);
    PolyMesh::Vertex *l20 = std::min(l2->v, l2->next->v);
    PolyMesh::Vertex *l21 = std::max(l2->v, l2->next->v);
    if(l10 == l20 && l11 == l21) return true;
    return false;
  }
};


void swap(double &a, double &b)
{
  double temp = a;
  a = b;
  b = temp;
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
      swap(xRed, xBlue);
      swap(yRed, yBlue);
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
      swap(xRed, xBlue);
      swap(yRed, yBlue);
      xBlue = -xBlue;
      yBlue = -yBlue;
      xRed = -xRed;
      yRed = -yRed;
      RESULT += 3 * BIG;
    }
    else
      printf("Hilbert failed %d %d", coordRed, coordBlue);
    BIG /= 4;
  }
  return RESULT;
}

static PolyMesh::Face *Walk(PolyMesh::Face *f, double x, double y)
{
  double POS[2] = {x, y};
  PolyMesh::HalfEdge *he = f->he;

  while(1) {
    PolyMesh::Vertex *v0 = he->v;
    PolyMesh::Vertex *v1 = he->next->v;
    PolyMesh::Vertex *v2 = he->next->next->v;

    double s0 = -robustPredicates::orient2d(v0->position, v1->position, POS);
    double s1 = -robustPredicates::orient2d(v1->position, v2->position, POS);
    double s2 = -robustPredicates::orient2d(v2->position, v0->position, POS);

    if(s0 >= 0 && s1 >= 0 && s2 >= 0) {
      /* printf("Face %g %g %g / %g %g %g / %g %g %g \n",
                v0->position.x(), v0->position.y(), v0->position.z(),
                v1->position.x(), v1->position.y(), v1->position.z(),
                v2->position.x(), v2->position.y(), v2->position.z());
                printf("point %g %g CURRENT FACE %p %g %g %g\n", x,y,he->f,
                s0,s1,s2); */
      // getchar();
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
                 "face %g %g %g / %g %g %g / %g %g %g "
                 "(orientation tests %g %g %g)", x, y,
                 v0->position.x(), v0->position.y(), v0->position.z(),
                 v1->position.x(), v1->position.y(), v1->position.z(),
                 v2->position.x(), v2->position.y(), v2->position.z(),
                 s0, s1, s2);
    }
    if(he == nullptr) break;
  }
  // should only come here wether the triangulated domain is not convex
  return nullptr;
}

static int delaunayEdgeCriterionPlaneIsotropic(PolyMesh::HalfEdge *he, void *)
{
  if(he->opposite == nullptr) return -1;
  PolyMesh::Vertex *v0 = he->v;
  PolyMesh::Vertex *v1 = he->next->v;
  PolyMesh::Vertex *v2 = he->next->next->v;
  PolyMesh::Vertex *v = he->opposite->next->next->v;

  // FIXME : should be oriented anyway !
  double result = -robustPredicates::incircle(v0->position, v1->position,
                                              v2->position, v->position);

  return (result > 0) ? 1 : 0;
}

extern "C" {
  int  polymesh_n_faces(const PolyMesh *pm);
  void polymesh_faces(const PolyMesh *pm, int *faces);
  PolyMesh *polymesh_new(double xmin[2], double xmax[2]);
  void polymesh_delete(PolyMesh *pm);
  void polymesh_add_points(PolyMesh *pm, int n, double *x, int *tags);
}

int polymesh_n_faces(const PolyMesh *pm) {
  return pm->faces.size();
}

void polymesh_faces(const PolyMesh *pm, int *faces) {
  for (size_t i = 0; i < pm->faces.size(); ++i) {
    PolyMesh::HalfEdge *he = pm->faces[i]->he;
    for (int j = 0; j < 3; ++j) {
      faces[i*3+j] = he->v->data;
      he = he->next;
    }
  }
}

void plymesh_delete(PolyMesh *pm) {
  delete pm;
}

static void get_bounding_box(int n, double *x, SVector3 &bbmin, SVector3 &bbmax) {
  bbmin = SVector3(std::numeric_limits<double>::max());
  bbmax = SVector3(-std::numeric_limits<double>::max());
  bbmin[2] = 0;
  bbmax[2] = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < 2; ++j) {
      bbmin[j] = std::min(x[i*2+j], bbmin[j]);
      bbmax[j] = std::max(x[i*2+j], bbmax[j]);
    }
  }
  for (int j = 0; j < 3; ++j) {
    double L = bbmax[j]-bbmin[j];
    bbmin[j] -= L*0.1;
    bbmax[j] += L*0.1;
  }
}

PolyMesh *polymesh_new(double xmin[2], double xmax[2]) {
  PolyMesh *pm = new PolyMesh;
  double d[2] = {xmax[0]-xmin[0],xmax[1]-xmin[1]};
  pm->initialize_rectangle(xmin[0]-d[0]*0.1, xmax[0]+d[0]*0.1,
      xmin[1]-d[1]*0.1, xmax[1]+d[1]*0.1);
  return pm;
}


void polymesh_add_points(PolyMesh *pm, int n, double *x, int *tags)
{
  std::vector<size_t> HC(n), IND(n);
  PolyMesh::Face *f = pm->faces[0];
  SVector3 bbmin(0), bbmax(0);
  get_bounding_box(n, x, bbmin, bbmax);
  SVector3 bbcenter = (bbmin+bbmax)/2;
  for(size_t i = 0; i < n; i++) {
    HC[i] = HilbertCoordinates(x[i*2], x[i*2+1], bbcenter.x(), bbcenter.y(),
                               bbmax.x() - bbcenter.x(), 0, 0,
                               bbmax.y() - bbcenter.y());
    IND[i] = i;
  }
  std::sort(IND.begin(), IND.end(),
            [&](size_t i, size_t j) { return HC[i] < HC[j]; });

  for(size_t i = 0; i < n; i++) {
    size_t I = IND[i];
    f = Walk(f, x[I*2], x[I*2+1]);
    pm->split_triangle(i, x[I*2], x[I*2+1], 0, f, delaunayEdgeCriterionPlaneIsotropic, nullptr);
    pm->vertices[pm->vertices.size() - 1]->data = tags[I];
  }
}
