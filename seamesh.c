#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef double (mesh_size_cb_t)(double x, double y,double z);

int find_neigh_id(int *t0, int *t1) {
  for (int i = 0; i<3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if(t0[i] == t1[j] && t0[(i+1)%3] == t1[(j+1)%3]) return i;
      if(t0[i] == t1[(j+1)%3] && t0[(i+1)%3] == t1[j]) return i;
    }
  }
  printf("neighbor not found\n");
  exit(1);
  return -1;
}

int is_inside(const double p[2], const double a[2], const double b[2], const double c[2]) {
  double dxdxi[2][2] = {{b[0]-a[0],c[0]-a[0]},{b[1]-a[1],c[1]-a[1]}};
  double det = dxdxi[0][0]*dxdxi[1][1]-dxdxi[1][0]*dxdxi[0][1];
  double dxidx[2][2] = {{dxdxi[1][1]/det,-dxdxi[0][1]/det},{-dxdxi[1][0]/det,dxdxi[1][1]/det}};
  double dp[2] = {p[0]-a[0],p[1]-a[1]};
  double xi[2] = {dxidx[0][0]*dp[0]+dxidx[0][1]*dp[1],dxidx[1][0]*dp[0]+dxidx[1][1]*dp[1]};
  return xi[0] >= -1e-8 && xi[1] >= -1e-8 && xi[0]+xi[1] <= 1+1e-8;
}

double triangle_signed_area(const double *a, const double *b, const double *c) {
  double dxdxi[2][2] = {{b[0]-a[0],c[0]-a[0]},{b[1]-a[1],c[1]-a[1]}};
  double det = dxdxi[0][0]*dxdxi[1][1]-dxdxi[1][0]*dxdxi[0][1];
  return det/2;
}

static double vdist(const double *v0, const double *v1) {
  return hypot(v1[0]-v0[0],v1[1]-v0[1]);
}

static double triangle_quality(const double *v0, const double *v1, const double *v2) {
  double a = vdist(v0,v1);
  double b = vdist(v0,v2);
  double c = vdist(v1,v2);
  double q=(a+b-c)*(a+c-b)*(b+c-a)/(2*a*b*c);
  if (q!=q) q = 0;
  return fabs(q);
}

static int getid3(const int t[3], int i) {
  if(t[0] == i) return 0;
  if(t[1] == i) return 1;
  if(t[2] == i) return 2;
  printf("element not found\n");
  return -1;
}


/*
 *          t3                    t3        
 *     i1 _______ i3         i1 _______ i3  
 *       |\_     |     ->      |     _/|    
 *   t4  |  \_t1 | t2  ->  t4  | t1_/  | t2 
 *       | t0 \_ |     ->      | _/ t0 |    
 *       |______\|     ->      |/______|    
 *     i2        i0          i2        i0   
 *          t5                    t5        
 */
static inline int swap(int t0, int e0, int *tris, double *x, int *neighbours) {
  int t1 = neighbours[t0*3+e0];
  int *tri0 = tris+t0*3;
  int *tri1 = tris+t1*3;
  int i0 = tri0[e0];
  int i1 = tri0[(e0+1)%3];
  int i2 = tri0[(e0+2)%3];
  int e1 = getid3(tris+t1*3,i1);
  int i3 = tri1[(e1+2)%3];
  if (triangle_signed_area(&x[i0*2],&x[i3*2],&x[i2*2])<=0) return 0;
  if (triangle_signed_area(&x[i1*2],&x[i2*2],&x[i3*2])<=0) return 0;
  double minqo = fmin(triangle_quality(&x[i0*2],&x[i1*2],&x[i2*2]),
                      triangle_quality(&x[i1*2],&x[i0*2],&x[i3*3]));
  double minqn = fmin(triangle_quality(&x[i0*2],&x[i3*2],&x[i2*2]),
                      triangle_quality(&x[i1*2],&x[i2*2],&x[i3*3]));
  if(minqn <= minqo*1.1) return 0;
  int t2 = neighbours[t1*3+(e1+1)%3];
  int t3 = neighbours[t1*3+(e1+2)%3];
  int t4 = neighbours[t0*3+(e0+1)%3];
  int t5 = neighbours[t0*3+(e0+2)%3];
  tri0[0] = i0;
  tri0[1] = i3;
  tri0[2] = i2;
  tri1[0] = i1;
  tri1[1] = i2;
  tri1[2] = i3;
  neighbours[t0*3+0] = t2;
  neighbours[t0*3+1] = t1;
  neighbours[t0*3+2] = t5;
  neighbours[t1*3+0] = t4;
  neighbours[t1*3+1] = t0;
  neighbours[t1*3+2] = t3;
  //i0->setF(t0);
  //i1->setF(t1);
  //modified[t0] = 1;
  //modified[t1] = 1;
  if (t2>=0) {
    neighbours[t2*3+getid3(neighbours+t2*3,t1)] = t0;
    //modified[t2] = 1;
  }
  if (t4>=0) {
    neighbours[t4*3+getid3(neighbours+t4*3,t0)] = t1;
    //modified[t4] = 1;
  }
  /*
  if (t3>=0) {
    t3->modified = 1;
  }*/
  return 1;
}

int *get_cavity(int t, int v, int *tri, int *neigh) {
  int cur = t;
  int n = 0;
  do {
    int iv = getid3(&tri[cur*3],v);
    n++;
    cur = neigh[cur*3+iv];
  }while (cur != -1 && cur != t);
  if (cur != t) {
    cur = t;
    int iv = getid3(&tri[cur*3],v);
    cur = neigh[cur*3+(iv+2)*3];
    while (cur != -1 && cur != t) {
      n++;
      int iv = getid3(&tri[cur*3],v);
      cur = neigh[cur*3+(iv+2)%3];
    }
  }
  int *r = malloc(sizeof(int)*n+1);
  n = 0;
  do {
    int iv = getid3(&tri[cur*3],v);
    r[n++] = cur;
    cur = neigh[cur*3+iv];
  }while (cur != -1 && cur != t);
  if (cur != t) {
    cur = t;
    int iv = getid3(&tri[cur*3],v);
    cur = neigh[cur*3+(iv+2)*3];
    while (cur != -1 && cur != t) {
      r[n++] = cur;
      int iv = getid3(&tri[cur*3],v);
      cur = neigh[cur*3+(iv+2)%3];
    }
  }
  r[n] = -1;
  return r;
}

static double optimal_boundary_collapse(int t0, int e0, int *tri, int *neigh, double *x, int **tri_color) {
  int v0 = tri[t0*3+e0];
  int v1 = tri[t0*3+(e0+1)%3];
  double d = vdist(&x[v0*2],&x[v1*2]);
  double d1=0, d0=0;
  int cur = t0;
  do {
    int iv1 = getid3(&tri[cur*3],v1);
    int n = neigh[cur*3+iv1];
    if (n == -1 || (*tri_color)[n] != (*tri_color)[cur]) {
      d1 = vdist(&x[v1*2],&x[tri[cur*3+(iv1+1)%3]*2]);
      break;
    }
    cur = n;;
  } while (cur != t0);
  cur = t0;
  do {
    int iv0 = getid3(&tri[cur*3],v0);
    int n = neigh[cur*3+(iv0+2)%3];
    if (n == -1 || (*tri_color)[n] != (*tri_color)[cur]) {
      d0 = vdist(&x[v0*2],&x[tri[cur*3+(iv0+2)%3]*2]);
      break;
    }
    cur = n;
  } while (cur != t0);
  double xi = fmin(1,fmax(0,(0.5*(d+d0+d1)-d0)/d));
  return 1-xi;
}

//substitute v0 by v1
void substitute_vertex(int t0, int v0, int v1, int *tri, int *neigh) {
  int cur = t0;
  do {
    int iv0 = getid3(&tri[cur*3],v0);
    if (cur != t0)
      tri[cur*3+iv0] = v1;
    cur = neigh[cur*3+iv0];
  } while (cur != t0 && cur != -1);
  if (cur == -1) {
    cur = t0;
    do {
      int iv0 = getid3(&tri[cur*3],v0);
      if (cur != t0)
        tri[cur*3+iv0] = v1;
      cur = neigh[cur*3+(iv0+2)%3];
    }while (cur != t0 && cur != -1);
  }
}

void substitute_neighbour(int t0, int t1, int rm, int *neigh) {
  if (t0 >= 0)
    neigh[t0*3+getid3(&neigh[t0*3],rm)] = t1;
  if (t1 >= 0)
    neigh[t1*3+getid3(&neigh[t1*3],rm)] = t0;
}

static int check_vertex_move(int t, int v, double *moveto, int ignore0, int ignore1, int *tri, int *neigh, double *x) {
  int cur = t;
  do {
    int iv = getid3(&tri[cur*3],v);
    if (cur != ignore0 && cur != ignore1) {
      if (triangle_signed_area(moveto,&x[tri[cur*3+(iv+1)%3]*2],&x[tri[cur*3+(iv+2)%3]*2])<0) return 0;
    }
    cur = neigh[cur*3+iv];
  } while (cur != t && cur != -1);
  if (cur == -1) {
    cur = t;
    do {
      int iv = getid3(&tri[cur*3],v);
      if (cur != ignore0 && cur != ignore1 && cur != t) {
        if (triangle_signed_area(moveto,&x[tri[cur*3+(iv+1)%3]*2],&x[tri[cur*3+(iv+2)%3]*2])<0) return 0;
      }
      cur = neigh[cur*3+(iv+2)%3];
    }while (cur != t && cur != -1); 
  }
  return 1;
}
/*
     v1                
     /|\
    / | \
   /  |  \             v0
  /   |   \           _/\_
 / t0 | t1 \        _/    \_
/_____|_____\      /        \
     v0
*/
// collapse tri[e0] on tri[(e0+1)%3] if possible.
static int collapse_edge(int t0, int e0, double xi, int *tri, int *neigh, double *x, double *length, int *vcolor) {
  int v0 = tri[t0*3+e0];
  int v1 = tri[t0*3+(e0+1)%3];
  double a = xi, b = 1-a;
  int t1 = neigh[t0*3+e0];
  double xmid[2] = {a*x[v0*2+0]+b*x[v1*2+0],a*x[v0*2+1]+b*x[v1*2+1]};
  if(!check_vertex_move(t0,v1,xmid,t0,t1,tri,neigh,x)) return 0;
  if(!check_vertex_move(t0,v0,xmid,t0,t1,tri,neigh,x)) return 0;
  x[v0*2+0] = xmid[0];
  x[v0*2+1] = xmid[1];
  length[v0] = a*length[v0] + b*length[v1];
  vcolor[v0] = vcolor[v0] > vcolor[v1] ? vcolor[v0] : vcolor[v1];
  //v0->setTag(std::max(v0->tag(),v1->tag()));
  substitute_vertex(t0,v1,v0,tri,neigh);
  substitute_neighbour(neigh[t0*3+(e0+1)%3],neigh[t0*3+(e0+2)%3],t0,neigh);
  tri[t0*3] = -1;
  if (t1 >= 0) {
    int e1 = getid3(&neigh[t1*3],t0);
    substitute_neighbour(neigh[t1*3+(e1+1)%3],neigh[t1*3+(e1+2)%3],t1,neigh);
    tri[t1*3] = -1;
  }
  return 1;
}

void gen_boundaries(int n_vertices, double *x, int *tag, int n_tri, int *tri, int *neigh_unsorted, double lon, double lat, mesh_size_cb_t *lc, double lcmin, int *n_tri_o, int **tri_o, int **tri_color, int *n_points_o_p, double **x_o_p)
{
  double *length = malloc(sizeof(double)*n_vertices);
  for (int i=0; i<n_vertices; ++i) {
    length[i] = lc(x[i*2+0],x[i*2+1],0);
  }
  //check orientation
  for (int i=0; i<n_tri; ++i){
    int *t = tri+i*3;
    if (triangle_signed_area(x+t[0]*2,x+t[1]*2,x+t[2]*2) < 0) {
      int sw = t[0];
      t[0] = t[1];
      t[1] = sw;
    }
  }
  //build neighbors
  int *neigh = malloc(sizeof(int)*n_tri*3);
  for (int i=0; i<n_tri; ++i){
    for (int j=0 ;j<3; ++j)
      neigh[i*3+j] = -1;
    for (int j = 0; j<3; ++j) {
      int n = neigh_unsorted[i*3+j];
      if (n<0) continue;
      neigh[i*3+find_neigh_id(tri+i*3,tri+n*3)] = n;
    }
  }
  // find vertex with first point
  int first = -1;
  double xp[2] = {lon,lat};
  for (int i=0; i<n_tri; ++i) {
    int *t = tri+i*3;
    if(is_inside(xp,x+t[0]*2,x+t[1]*2,x+t[2]*2)) {
      first = i;
    }
  }
  if(first==-1) {
    printf("no triangle contains the first point\n");
    exit(1);
  }

  // color triangles
  *tri_color = malloc(sizeof(int)*n_tri);
  for (int i=0; i<n_tri; ++i) {
    (*tri_color)[i] = 0;
  }
  int *stack = malloc(sizeof(int)*n_tri);
  int stack_p = 0;
  stack[0] = first;
  (*tri_color)[first] = 1;
  while(stack_p >= 0) {
    int c = stack[stack_p--];
    for (int j=0; j<3; ++j) {
      int n = neigh[c*3+j];
      if (n >=0 && !(*tri_color)[n]) {
        int i0 = tri[c*3+j];
        int i1 = tri[c*3+(j+1)%3];
        double *p0 = x+2*i0;
        double *p1 = x+2*i1;
        double d2 = (p1[0]-p0[0])*(p1[0]-p0[0])+(p1[1]-p0[1])*(p1[1]-p0[1]);
        double lmin = fmin(length[i0],length[i1]);
        if(d2 > lmin*lmin) {
          (*tri_color)[n] = 1;
          stack[++stack_p] = n;
        }
      }
    }
  }
  
  free(stack);
  int *vcolor = malloc(sizeof(int)*n_vertices);
  for (int i=0; i<n_vertices; ++i) {
    vcolor[i] = 0;
  }
  for (int i=0; i<n_tri; ++i) {
    if ((*tri_color)[i] == 1) {
      for (int j=0; j<3; ++j) {
        vcolor[tri[i*3+j]] = 1;
      }
    }
  }

  //swap - collapses
  size_t nops ;
  int iter = 0;
  do {
    iter += 1;
    nops = 0;
    for (int i=0;i<n_tri;i++) {
      int *t = tri+i*3;
      for (int j=0; j<3 && t[0] != -1;++j) {
        int n = neigh[i*3+j];
        if (n != -1 && (*tri_color)[n] == (*tri_color)[i]){
          //nops += swap(i,j,tri,x,neigh);
        }
        int i0 = t[j];
        int i1 = t[(j+1)%3];
        if((*tri_color)[i] == 0 && (n == -1 || (*tri_color)[n] !=1)) {
          // we are outside the domain
          // and so is our neighbour (if any)
          // remove point that does not touch the boundary
          if (vcolor[i1] == 0) {
            int r = collapse_edge(i,j,1,tri,neigh,x,length,vcolor);
            nops += r;
          }
          else if (vcolor[i0] == 0){
            int r = collapse_edge(i,j,0,tri,neigh,x,length,vcolor);
            nops += r;
          }
        }
        if((*tri_color)[i] == 1 && (n == -1 || (*tri_color)[n] ==0)
            && vdist(&x[i0*2],&x[i1*2]) < 0.7*fmin(length[i0],length[i1])) {
          // we are on the boundary
          // collapse edge if too short (and possible)
              nops += 
                collapse_edge(i,j,optimal_boundary_collapse(i,j,tri,neigh,x,tri_color),tri,neigh,x,length,vcolor) ||
                collapse_edge(i,j,0,tri,neigh,x,length,vcolor) || collapse_edge(i,j,1,tri,neigh,x,length,vcolor);
        }
      }
    }
    printf("nops = %i\n",nops);
  } while (nops != 0 && iter < 100);
  for (int i = 0; i < n_tri; ++i) {
    int *t = &tri[i*3];
          if(t[0] == -1) continue;
    if(triangle_signed_area(&x[t[0]*2],&x[t[1]*2],&x[t[2]*2])<0) {
      printf("invalid triangle !!\n");
    }
  }

  free(neigh);
  free(length);
  free(vcolor);
  int *keep_vertices = malloc(sizeof(int)*n_vertices);
  for (int i = 0; i < n_vertices; ++i)
    keep_vertices[i] = -1;
  int n_vertices_o = 0;
  *n_tri_o = 0;
  for (int i = 0; i < n_tri; ++i) {
    int *t = tri+i*3;
    if (t[0] == -1) continue;
    *n_tri_o += 1;
    for (int j = 0; j < 3; ++j){
      if (keep_vertices[t[j]] == -1) {
        n_vertices_o += 1;
        keep_vertices[t[j]] = 1;
      }
    }
  }
  double *x_o = malloc(sizeof(double)*2*n_vertices_o);
  n_vertices_o = 0;
  for (int i = 0; i < n_vertices; ++i) {
    if(keep_vertices[i] == 1) {
      keep_vertices[i] = n_vertices_o;
      x_o[n_vertices_o*2+0] = x[i*2+0];
      x_o[n_vertices_o*2+1] = x[i*2+1];
     n_vertices_o += 1;
    }
  }
  *tri_o = malloc((*n_tri_o)*sizeof(int)*3);
  int *tri_color_o = malloc(sizeof(int)*(*n_tri_o));
  int i_tri_o = 0;
  for (int i = 0; i < n_tri; ++i) {
    int *t = tri+i*3;
    if (t[0] == -1) continue;
    int *t_o = *tri_o+i_tri_o*3;
    for (int j = 0; j<3; ++j) {
      t_o[j] = keep_vertices[t[j]];
    }
    tri_color_o[i_tri_o] = (*tri_color)[i];
    i_tri_o += 1;
  }
  free(*tri_color);
  *tri_color = tri_color_o;
  free(keep_vertices);
  *n_points_o_p = n_vertices_o;
  *x_o_p = x_o;
  printf("end gen boundaries\n");
}

static int hecmp(const void *p0, const void *p1) {
  const int *i0 = (const int *)p0;
  const int *i1 = (const int *)p1;
  return i0[0] == i1[0] ? i0[1]-i1[1] : i0[0]-i1[0];
}

int *build_neighbours(int n_tri, const int *tri) {
  int *hedges = malloc(sizeof(int)*4*n_tri*3);
  for (int i = 0; i < n_tri; ++i) {
    const int *t = tri+i*3;
    for (int j = 0; j < 3; ++j) {
      int ihe = i*3+j;
      int *he = hedges + ihe*4;
      int i1 = (i+1)%3;
      if(t[i1] < t[i]) {
        he[0] = t[i1];
        he[1] = t[i];
      }
      else {
        he[0] = t[i];
        he[1] = t[i1];
      }
      he[2] = i;
      he[3] = j;
    }
  }
  qsort(hedges,n_tri*3,4*sizeof(int),hecmp);
  int *neighbours = malloc(sizeof(int)*n_tri*3);
  for (int i = 0; i < n_tri*3; ++i) {
    neighbours[i] = -1;
  }
  for (int i = 0; i+1<n_tri*3; ++i) {
    int *he0 = hedges+i*4;
    int *he1 = hedges+i*4+4;
    if (hecmp(he0,he1) == 0) {
      neighbours[he0[2]*3+he0[3]] = he1[2];
      neighbours[he1[2]*3+he1[3]] = he0[2];
    }
  }
  free(hedges);
  return neighbours;
}

void write_geo(const char *filename, double lc, int n_tri, const int *tri, const int *tri_color, int n_points, const double *x_p, const int *vtag/*, std::map<int,std::string> boundary_names*/)
{
  int *neighbours = build_neighbours(n_tri, tri);
  FILE *f = fopen(filename,"w");
  //std::map<std::string,std::vector<int> > physicals;
  fprintf(f,"mesh.lcintegrationprecision = 1e-4;\n");
  fprintf(f,"lc=%g;\n",lc);
  size_t ip = 1;
  size_t il = 1;
  size_t ill = 1;
  char *touched = malloc(sizeof(char)*n_tri);
  for (int i = 0; i < n_tri; ++i) {
    touched[i] = 0;
  }
  for (int i=0;i<n_tri;i++){
    const int *t = tri+3*i;
    if(tri_color[i] != 0 || touched[i]) continue;
    for (int j=0; j<3; ++j) {
      int neigh = neighbours[i*3+j];
      if (neigh == -1 || tri_color[neigh] || touched[neigh]) continue;
      int firstp = ip;
      int firstpll = ip;
      int firstl = il;
      int firstv = t[j], curv = t[j];
      int curt = i;
      int iv = j;
      int oldtag = -1;
      do {
        touched[curt] = 1;
        fprintf(f,"point(%i)={%.16g,%.16g,%.16g};\n",ip++,x_p[curv*3+0],x_p[curv*3+1],x_p[curv*3+2]);
        int tag = vtag[curv];
        if (oldtag != -1 && oldtag != tag) {
          //physicals[boundary_names[oldtag]].push_back(il);
          fprintf(f,"spline(%i)={%i:%i};\n",il++,firstp,ip-1);
          firstp = ip-1;
        }
        oldtag = tag;
        iv = (iv+1)%3;
        curv = t[iv];
        while(tri_color[neighbours[curt*3+iv]] == 0) {
          curt = neighbours[curt*3+iv];
          iv = getid3(tri+curt*3,curv);
        }
      }while(curv != firstv);
      //physicals[boundary_names[oldtag]].push_back(il);
      fprintf(f,"spline(%i)={%i:%i,%i};\n",il++,firstp,ip-1,firstpll);
      fprintf(f,"line loop(%i={%i:%i};\n",firstl,il-1);
      il++;
    }
  }
  /*f << "physical surface(\"domain\") = {1};\n";
  for(auto p: physicals) {
    f << "physical line(\"" << p.first << "\") = {";
    for(size_t i = 0; i < p.second.size()-1; ++i) {
      f << p.second[i] << ",";
    }
    f << p.second.back() << "};\n";

  }*/
  fprintf(f,"Plane Surface(1)={1:%i};\n",ill-1);
  fprintf(f,"Mesh.CharacteristicLengthExtendFromBoundary = 0;\n");
  fprintf(f,"Mesh.CharacteristicLengthFromPoints = 0;\n");
  fclose(f);
  free(neighbours);
}
