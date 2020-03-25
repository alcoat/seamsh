import numpy as np
import ctypes as c
from ctypes.util import find_library
import os
import platform 
from gdal import ogr,osr

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

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

def read_points_from_shp(filename,lcmin,ellps) :
    points = []
    tags = []
    proj = osr.SpatialReference()
    #proj.ImportFromProj4("+ellps=WGS84 +proj=utm +utm_zone=31")
    proj.ImportFromProj4("+ellps=WGS84 +proj=geocent")
    driver = ogr.GetDriverByName('ESRI Shapefile')
    data = driver.Open(filename,0)
    layer = data.GetLayer()
    layerproj = layer.GetSpatialRef()
    reprojector = osr.CoordinateTransformation(layerproj,proj)
    for i in layer :
        try :
            phys = i.GetField("entity")
        except :
            phys = -1
        g = i.GetGeometryRef()
        lastp = None
        pts = np.array(reprojector.TransformPoints(g.GetPoints()))
        d=pts[1:,:]-pts[:-1,:]
        l=np.linalg.norm(d,axis=1)
        n=np.ceil(l/lcmin).astype(np.int32)
        for i in range(pts.shape[0]-1) :
            for j in range(n[i]):
                alpha = j/n[i]
                points.append((1-alpha)*pts[i]+alpha*pts[i+1])
                tags.append(phys)
        points.append(pts[-1])
        tags.append(phys)
    return np.array(points),np.array(tags)


def np2c(a,dtype=np.float64) :
    tmp = np.require(a,dtype,"C")
    r = c.c_void_p(tmp.ctypes.data)
    r.tmp = tmp
    return r

def gen_boundaries(x,tag,tri,lon,lat,mesh_size_cb) :
    n_ll = c.c_int()
    n_l = c.c_int()
    n_xo = c.c_int()
    p_xo = c.POINTER(c.c_double)()
    p_l = c.POINTER(c.c_int)()
    p_ll = c.POINTER(c.c_int)()
    mesh_size = np.array([mesh_size_cb(xi[0],xi[1],0)for xi in x])
    cb = c.CFUNCTYPE(c.c_double,c.c_double,c.c_double,c.c_double)(mesh_size_cb)
    lib.gen_boundaries_from_points(
            c.c_int(x.shape[0]),np2c(x),np2c(tag,np.int32),
            c.c_int(tri.shape[0]),np2c(tri,np.int32),
            c.c_double(lon),c.c_double(lat),
            np2c(mesh_size),
            c.byref(p_xo),c.byref(n_xo),
            c.byref(p_l),c.byref(n_l),
            c.byref(p_ll),c.byref(n_ll))
    xo = np.ctypeslib.frombuffer(c.cast(p_xo,c.POINTER(n_xo.value*2*c.c_double)).contents,dtype=np.float64).reshape([-1,2]).copy()
    lib.libcfree(p_xo)
    l = np.ctypeslib.frombuffer(c.cast(p_l,c.POINTER(n_l.value*3*c.c_int)).contents,dtype=np.int32).reshape([-1,3]).copy()
    lib.libcfree(p_l)
    ll = np.ctypeslib.frombuffer(c.cast(p_ll,c.POINTER(n_ll.value*2*c.c_int)).contents,dtype=np.int32).reshape([-1,2]).copy()
    lib.libcfree(p_ll)
    return xo,l,ll
