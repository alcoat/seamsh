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

static void merge_tag_list(int *tag, int v0, int v1) {
  // attach v1 chain at the end of v0 chain
  int v0endtag = v0;
  while(tag[v0endtag*2+1] != -1) v0endtag = tag[v0endtag*2+1];
  tag[v0endtag*2+1] = v1;
  return;
  // remove duplicates
  int end = v0;
  int last = v0;
  int cur;
  for (cur=v0; cur != -1; cur = tag[cur*2+1]) {
    int found = 0;
    for (int i = v0; i!= end; i = tag[i*2+1]) {
      found |= tag[i*2] == tag[cur*2];
    }
    if (!found) {
      tag[end*2] = tag[cur*2];
      last = end;
      end = tag[end*2+1];
    }
  }
  if (end!=-1) {
    tag[last*2+1] = -1;
  }
}

static int has_tag(const int *tag, int v, int t) {
  while (v != -1) {
    if (tag[v*2] == t) return 1;
    v = tag[v*2+1];
  }
  return 0;
}

static int get_segment_tag(const int *tag, int v0, int v1, int current_tag) {
  if (has_tag(tag,v0,current_tag) && has_tag(tag,v1,current_tag))
    return current_tag;
  for (int i = v0; i!= -1; i = tag[i*2+1]) {
    if(has_tag(tag,v1,tag[i*2])) return tag[i*2];
  }
  if (has_tag(tag,v0,current_tag) || has_tag(tag,v1,current_tag))
    return current_tag;
  return tag[v0*2];
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
static int collapse_edge(int t0, int e0, double xi, int *tri, int *neigh, double *x, double *length, int *vcolor,int *tag) {
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
  merge_tag_list(tag,v0,v1);
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

void *vector_append(void **p, int *psize, int elem_size) {
  if (*psize == 0) {
    *p = malloc(elem_size*2);
  }
  else if (*psize%2==0) {
    *p = realloc(*p,(*psize)*2*elem_size);
  }
  *psize += 1;
  return ((char*)*p) + ((*psize-1)*elem_size);
}

static int coord_vector_append(double **pts, int *npts, const double *x) {
  double *xn = (double*)vector_append((void**)pts,npts,2*sizeof(double));
  xn[0] = x[0];
  xn[1] = x[1];
  return (*npts)-1;
}

static int int_vector_append(int **v, int *n, int i) {
  return *(int*)vector_append((void**)v,n,sizeof(int)) = i;
}

static int double_vector_append(double **v, int *n, double i) {
  return *(double*)vector_append((void**)v,n,sizeof(double)) = i;
}

static void extract_boundaries(const double *x_p, const double *x_size, const int *vtag, int n_tri, const int *tri, const int *tri_color,double **xo, int *nxo, int **l, int *nl)
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
  int *touched = malloc(sizeof(int)*n_tri_f);
  for(int i = 0; i < n_tri_f; ++i)
    touched[i] = 0;
  // build line loops
  int *ll = NULL;
  int nll = 0;

  for(int i = 0; i < n_tri_f; ++i){
    if (touched[i] == 1) continue;
    for (int j = 0; j<3; ++j) {
      if(neighbours[i*3+j]>=0) continue;
      int firstp = tri_f[i*3+j];
      int curt = i;
      int p = j;
      touched[i] = 1;
      int curp;
      do {
        curp = tri_f[curt*3+(p+1)%3];
        int_vector_append(&ll,&nll,curp);
        p = getid3(tri_f+curt*3,curp);
        while (neighbours[curt*3+p] >=0) {
          curt = neighbours[curt*3+p]; 
          p = getid3(tri_f+curt*3,curp);
        }
        touched[curt] = 1;
      } while(curp != firstp || curt !=i);
      int_vector_append(&ll,&nll,-1);
      break;
    }
  }
  // build lines
  *xo = NULL;
  *l = NULL;
  *nxo = 0, *nl = 0;
  double *osize = NULL;
  int nosize = 0;
  for (int i = 0; i < nll; ++i) {
    int *lli = ll+i;
    int length = 0;
    for (length = 0; lli[length] != -1; ++length);
    int start = 0;
    int ctag = get_segment_tag(vtag,lli[0],lli[1],-1);
    for (int j=1; j < length; ++j) {
      int ntag = get_segment_tag(vtag,lli[j],lli[(j+1)%length],ctag);
      if (ctag != ntag) {
        start = j;
        break;
      }
    }
    int firstp = *nxo;;
    ctag = get_segment_tag(vtag,lli[start],lli[(start+1)%length],-1);
    int_vector_append(l,nl,ctag);
    for (int j = 0; j < length; ++j) {
      int p = lli[(start+j)%length];
      int pnext = lli[(start+j+1)%length];
      int stag = get_segment_tag(vtag,p,pnext,ctag);
      coord_vector_append(xo,nxo,x_p+p*2);
      double_vector_append(&osize,&nosize,x_size[p]);
      int_vector_append(l,nl,(*nxo)-1);
      if (stag != ctag) {
        int_vector_append(l,nl,-1);
        int_vector_append(l,nl,stag);
        int_vector_append(l,nl,(*nxo)-1);
        ctag = stag;
      }
    }
    int_vector_append(l,nl,firstp);
    int_vector_append(l,nl,-1);
    i += length;
    // slightly extrude the boundary towards the domain
    double *normals = malloc(sizeof(double)*2*length);
    for (int j = 0; j < length; ++j) {
      const double *x0 = (*xo) + (firstp+j)*2;
      const double *x1 = (*xo) + (firstp+(j+length-1)%length)*2;
      double *n = normals + j*2;
      n[0] = x1[1] - x0[1];
      n[1] = -x1[0] + x0[0];
      double l = hypot(n[0], n[1]);
      n[0] /= l;
      n[1] /= l;
    }
    for (int j = 0; j < length; ++j) {
      const double *n0 = &normals[j*2];
      const double *n1 = &normals[((j+length-1)%length)*2];
      (*xo)[(firstp+j)*2+0] += osize[firstp+j]*(n0[0]+n1[0])/(50);
      (*xo)[(firstp+j)*2+1] += osize[firstp+j]*(n0[1]+n1[1])/(50);
    }
    free(normals);
  }
  free(osize);
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

static int *color_triangles(const double *x, int n_tri, const int *tri, const int *neigh, int first, const double *length) {
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

static void optimize_mesh(double *x, double *length, int nx, int n_tri, int *tri, int *tri_color, int *neigh,int *tag) {
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
            int r = collapse_edge(i,j,1,tri,neigh,x,length,vcolor,tag);
            nops += r;
          }
          else if (vcolor[i0] == 0){
            int r = collapse_edge(i,j,0,tri,neigh,x,length,vcolor,tag);
            nops += r;
          }
        }
        if(tri_color[i] == 1 && (n == -1 || tri_color[n] ==0)
            && vdist(&x[i0*2],&x[i1*2]) < 0.7*fmin(length[i0],length[i1])) {
          // we are on the boundary
          // collapse edge if too short (and possible)
          double opt = optimal_boundary_collapse(i,j,tri,neigh,x,tri_color);
          nops += 
            collapse_edge(i,j,opt,tri,neigh,x,length,vcolor,tag)
            || collapse_edge(i,j,0,tri,neigh,x,length,vcolor,tag)
            || collapse_edge(i,j,1,tri,neigh,x,length,vcolor,tag);
        }
      }
    }
  } while (nops != 0 && iter < 100);
  free(vcolor);
}

void gen_boundaries_from_points(int n_vertices, double *x, int *tag, int n_tri, int *tri, int first, double *mesh_size, double **xo, int *nxo, int **l, int *nl)
{
  orient_triangles(x,n_tri,tri);
  int *neigh = build_neighbours(n_tri,tri);
  int *tri_color = color_triangles(x,n_tri,tri,neigh,first,mesh_size);
  optimize_mesh(x, mesh_size, n_vertices, n_tri, tri, tri_color, neigh,tag);
  extract_boundaries(x, mesh_size, tag, n_tri, tri, tri_color, xo,  nxo, l, nl);
  free(tri_color);
  free(neigh);
}

void libcfree(void *p) {
  free(p);
}
