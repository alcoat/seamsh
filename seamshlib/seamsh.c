/*
* seamsh - Copyright (C) <2010-2020>
* <Universite catholique de Louvain (UCL), Belgium
* 	
* List of the contributors to the development of seamsh: see AUTHORS file.
* Description and complete License: see LICENSE file.
* 	
* This program (seamsh) is free software: 
* you can redistribute it and/or modify it under the terms of the GNU Lesser General 
* Public License as published by the Free Software Foundation, either version
* 3 of the License, or (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with this program (see COPYING file).  If not, 
* see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void *vector_push(void **p, int *psize, int elem_size) {
  if (*psize == 0) {
    *p = malloc(elem_size*2);
  }
  else if (*psize%2==0) {
    *p = realloc(*p,(*psize)*2*elem_size);
  }
  *psize += 1;
  return ((char*)*p) + ((*psize-1)*elem_size);
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
  exit(1);
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

static double optimal_boundary_collapse(int t0, int e0, int *tri, int *neigh, double *x, int *tri_color) {
  int v0 = tri[t0*3+e0];
  int v1 = tri[t0*3+(e0+1)%3];
  double d = vdist(&x[v0*2],&x[v1*2]);
  double d1=0, d0=0;
  int cur = t0;
  do {
    int iv1 = getid3(&tri[cur*3],v1);
    int n = neigh[cur*3+iv1];
    if (n == -1 || tri_color[n] != tri_color[cur]) {
      d1 = vdist(&x[v1*2],&x[tri[cur*3+(iv1+1)%3]*2]);
      break;
    }
    cur = n;;
  } while (cur != t0);
  cur = t0;
  do {
    int iv0 = getid3(&tri[cur*3],v0);
    int n = neigh[cur*3+(iv0+2)%3];
    if (n == -1 || tri_color[n] != tri_color[cur]) {
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

static int hecmp(const void *p0, const void *p1) {
  const int *i0 = (const int *)p0;
  const int *i1 = (const int *)p1;
  return i0[0] == i1[0] ? i0[1]-i1[1] : i0[0]-i1[0];
}

static int *build_neighbours(int n_tri, const int *tri) {
  int *hedges = malloc(sizeof(int)*4*n_tri*3);
  int ihe =0;
  for (int i = 0; i < n_tri; ++i) {
    const int *t = tri+i*3;
    if (t[0]<0) continue;
    for (int j = 0; j < 3; ++j) {
      int *he = hedges + ihe*4;
      int j1 = (j+1)%3;
      if(t[j1] < t[j]) {
        he[0] = t[j1];
        he[1] = t[j];
      }
      else {
        he[0] = t[j];
        he[1] = t[j1];
      }
      he[2] = i;
      he[3] = j;
      ihe++;
    }
  }
  qsort(hedges,ihe,4*sizeof(int),hecmp);
  int *neighbours = malloc(sizeof(int)*n_tri*3);
  for (int i = 0; i < n_tri*3; ++i) {
    neighbours[i] = -1;
  }
  for (int i = 0; i+1<ihe; ++i) {
    int *he0 = hedges+i*4;
    int *he1 = hedges+i*4+4;
    if (hecmp(he0,he1) == 0) {
      neighbours[he0[2]*3+he0[3]] = he1[2];
      neighbours[he1[2]*3+he1[3]] = he0[2];
      i++;
    }
  }
  free(hedges);
  return neighbours;
}

static void extract_boundaries(const double *x_p, const int *vtag, int n_tri, const int *tri, const int *tri_color, /*, std::map<int,std::string> boundary_names*/double **xo, int *nxo, int **l, int *nl, int **ll, int *nll)
{
  // filter triangles
  int n_tri_f = 0;
  int *tri_f =  malloc(sizeof(int)*n_tri*3);
  for (int i= 0; i < n_tri; ++i) {
    if (tri[i*3] < 0 || tri_color[i] != 1) continue;
    for (int j = 0; j < 3; ++j) {
      tri_f[n_tri_f*3+j] = tri[i*3+j];
    }
    n_tri_f++;
  }
  int *neighbours = build_neighbours(n_tri_f,tri_f);
  // walk
  *xo = NULL;
  *ll = NULL;
  *l = NULL;
  *nxo = 0, *nll = 0, *nl = 0;
  int *touched = malloc(sizeof(int)*n_tri_f);
  for(int i = 0; i < n_tri_f; ++i)
    touched[i] = 0;
  for(int i = 0; i < n_tri_f; ++i){
    if (touched[i] == 1) continue;
    for (int j = 0; j<3; ++j) {
      if(neighbours[i*3+j]>=0) continue;
      int firstp = tri_f[i*3+j];
      double *X = (double*)vector_push((void**)xo,nxo,sizeof(double)*2);
      X[0] = x_p[firstp*2+0];
      X[1] = x_p[firstp*2+1];
      int curp = tri_f[i*3+(j+1)%3];
      int curt = i;
      int *line = (int*)vector_push((void**)l,nl,3*sizeof(int));
      line[0] = (*nxo)-1;
      touched[i] = 1;
      while (curp != firstp) {
        double *X = (double*)vector_push((void**)xo,nxo,sizeof(double)*2);
        X[0] = x_p[curp*2+0];
        X[1] = x_p[curp*2+1];
        int p = getid3(tri_f+curt*3,curp);
        while (neighbours[curt*3+p] >=0) {
          curt = neighbours[curt*3+p]; 
          p = getid3(tri_f+curt*3,curp);
        }
        touched[curt] = 1;
        curp = tri_f[curt*3+(p+1)%3];
      }
      line[1] = (*nxo)-1;
      line[2] = line[0];
      break;
    }
  }
  free(neighbours);
  free(touched);
  free(tri_f);
}

static void orient_triangles(const double *x, int n_tri, int *tri) {
  for (int i=0; i<n_tri; ++i){
    int *t = tri+i*3;
    if (triangle_signed_area(x+t[0]*2,x+t[1]*2,x+t[2]*2) < 0) {
      int sw = t[0];
      t[0] = t[1];
      t[1] = sw;
    }
  }
}

static int find_triangle_containing_point(const double *x, int n_tri, const int *tri, double lon, double lat) {
  // find triangle with first point
  int first = -1;
  double xp[2] = {lon,lat};
  for (int i=0; i<n_tri; ++i) {
    const int *t = tri+i*3;
    if(is_inside(xp,x+t[0]*2,x+t[1]*2,x+t[2]*2)) {
      first = i;
    }
  }
  if(first==-1) {
    printf("no triangle contains the first point %g %g\n",lon,lat);
    exit(1);
  }
  return first;
}

static int *color_triangles(const double *x, int n_tri, const int *tri, const int *neigh, double lon, double lat, const double *length) {
  int first = find_triangle_containing_point(x,n_tri,tri,lon,lat);
  int *tri_color = malloc(sizeof(int)*n_tri);
  for (int i=0; i<n_tri; ++i) {
    tri_color[i] = 0;
  }
  int *stack = malloc(sizeof(int)*n_tri);
  int stack_p = 0;
  stack[0] = first;
  tri_color[first] = 1;
  while(stack_p >= 0) {
    int c = stack[stack_p--];
    for (int j=0; j<3; ++j) {
      int n = neigh[c*3+j];
      if (n >=0 && !tri_color[n]) {
        int i0 = tri[c*3+j];
        int i1 = tri[c*3+(j+1)%3];
        const double *p0 = x+2*i0;
        const double *p1 = x+2*i1;
        double d2 = (p1[0]-p0[0])*(p1[0]-p0[0])+(p1[1]-p0[1])*(p1[1]-p0[1]);
        double lmin = fmin(length[i0],length[i1]);
        if(d2 > lmin*lmin) {
          tri_color[n] = 1;
          stack[++stack_p] = n;
        }
      }
    }
  }
  free(stack);
  return tri_color;
}

static void optimize_mesh(double *x, double *length, int nx, int n_tri, int *tri, int *tri_color, int *neigh) {
  int *vcolor = malloc(sizeof(int)*nx);
  for (int i=0; i<nx; ++i) {
    vcolor[i] = 0;
  }
  for (int i=0; i<n_tri; ++i) {
    if (tri_color[i] == 1) {
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
        if (n != -1 && tri_color[n] == tri_color[i]){
          //nops += swap(i,j,tri,x,neigh);
        }
        int i0 = t[j];
        int i1 = t[(j+1)%3];
        if(tri_color[i] == 0 && (n == -1 || tri_color[n] !=1)) {
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
        if(tri_color[i] == 1 && (n == -1 || tri_color[n] ==0)
            && vdist(&x[i0*2],&x[i1*2]) < 0.7*fmin(length[i0],length[i1])) {
          // we are on the boundary
          // collapse edge if too short (and possible)
          nops += 
            collapse_edge(i,j,optimal_boundary_collapse(i,j,tri,neigh,x,tri_color),tri,neigh,x,length,vcolor) ||
            collapse_edge(i,j,0,tri,neigh,x,length,vcolor) || collapse_edge(i,j,1,tri,neigh,x,length,vcolor);
        }
      }
    }
  } while (nops != 0 && iter < 100);
  free(vcolor);
}

void gen_boundaries_from_points(int n_vertices, double *x, int *tag, int n_tri, int *tri, double lon, double lat, double *mesh_size, double **xo, int *nxo, int **l, int *nl, int **ll, int *nll)
{
  orient_triangles(x,n_tri,tri);
  int *neigh = build_neighbours(n_tri,tri);
  int *tri_color = color_triangles(x,n_tri,tri,neigh,lon,lat,mesh_size);
  optimize_mesh(x, mesh_size, n_vertices, n_tri, tri, tri_color, neigh);
  extract_boundaries(x, tag, n_tri, tri, tri_color/*, std::map<int,std::string> boundary_names*/,xo,  nxo, l,nl,ll,nll);
  free(tri_color);
  free(neigh);
}

void libcfree(void *p) {
  free(p);
}
