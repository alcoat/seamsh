import numpy as np
import ctypes as c
from ctypes.util import find_library
import os
import platform 
libdir = os.path.dirname(os.path.realpath(__file__))
if platform.system() == "Windows":
    libpath = os.path.join(libdir, "unref.dll")
elif platform.system() == "Darwin":
    libpath = os.path.join(libdir, "libunref.dylib")
else:
    libpath = os.path.join(libdir, "build/libseamesh.so")

lib = c.CDLL(libpath)

import numpy as np
import ctypes as c

def read_shp(filename, lcmin, ellps) :
    n = c.c_int()
    xp = c.POINTER(c.c_double)()
    tagp = c.POINTER(c.c_int)()
    lib.read_shp(filename.encode(),c.c_double(lcmin),ellps.encode(),c.byref(n),c.byref(xp),c.byref(tagp))
    x = np.ctypeslib.frombuffer(c.cast(xp,c.POINTER(n.value*3*c.c_double)).contents).reshape([-1,3]).copy()
    tag = np.ctypeslib.frombuffer(c.cast(tagp,c.POINTER(n.value*c.c_int)).contents,dtype=np.int32).copy()
    lib.libcfree(xp)
    lib.libcfree(tagp)
    return x,tag

def np2c(a,dtype=np.float64) :
    tmp = np.require(a,dtype,"C")
    r = c.c_void_p(tmp.ctypes.data)
    r.tmp = tmp
    return r

def gen_boundaries(x,tag,tri,lon,lat,mesh_size_cb,lcmin) :
    n_tri_p = c.c_int()
    n_points_p = c.c_int()
    x_p = c.POINTER(c.c_double)()
    tri_o_p = c.POINTER(c.c_int)()
    tri_color_p = c.POINTER(c.c_int)()
    cb = c.CFUNCTYPE(c.c_double,c.c_double,c.c_double,c.c_double)(mesh_size_cb)
    lib.gen_boundaries(
            c.c_int(x.shape[0]),np2c(x),np2c(tag,np.int32),
            c.c_int(tri.simplices.shape[0]),np2c(tri.simplices,np.int32),np2c(tri.neighbors,np.int32),
            c.c_double(lon),c.c_double(lat),
            cb,c.c_double(lcmin),
            c.byref(n_tri_p),c.byref(tri_o_p),c.byref(tri_color_p),
            c.byref(n_points_p), c.byref(x_p))
    tri_o = np.ctypeslib.frombuffer(c.cast(tri_o_p,c.POINTER(n_tri_p.value*3*c.c_int)).contents,dtype=np.int32).reshape([-1,3]).copy()
    tri_color_o = np.ctypeslib.frombuffer(c.cast(tri_color_p,c.POINTER(n_tri_p.value*c.c_int)).contents,dtype=np.int32).copy()
    x_o = np.ctypeslib.frombuffer(c.cast(x_p,c.POINTER(n_points_p.value*2*c.c_double)).contents,dtype=np.float64).reshape([-1,2]).copy()
    lib.libcfree(tri_o_p)
    lib.libcfree(tri_color_p)
    lib.libcfree(x_p)
    return x_o,tri_o,tri_color_o

def write_geo(filename,lc, x, tag, tri, tri_color) :
    lib.write_geo(filename.encode(),c.c_double(lc),np2c(x),np2c(tag,np.int32),c.c_int(tri.shape[0]),np2c(tri,np.int32),np2c(tri_color,np.int32))
#void write_geo(const char *filename, double lc, int n_tri, const int *tri, const int *tri_color, int n_points, const double *x_p, const int *vtag/*, std::map<int,std::string> boundary_names*/)
