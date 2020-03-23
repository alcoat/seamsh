#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "Vertex.h"
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include "proj_api.h"
#define ERROR std::cout << "wrong file " << __LINE__ << std::endl
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846	
#endif
static void swapIntByte(int &x) {
   x=((x & 0xff000000) >> 24) |
     ((x & 0x00ff0000) >> 8 ) |
     ((x & 0x0000ff00) << 8 ) |
     ((x & 0x000000ff) << 24); 
}

static void radiusScale(double R, double xyz[3]) {
  double f = R/sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  xyz[0] *= f;
  xyz[1] *= f;
  xyz[2] *= f;
}

void saturateEdge(std::vector<Vertex*> &vertices, const Vertex &v0, const Vertex &v1, double lcmin, int tag, bool withShiftedStart=false, bool withShiftedEnd=false)
{
  const double dx = v1.x()-v0.x();
  const double dy = v1.y()-v0.y();
  const double dz = v1.z()-v0.z();
  const double r0 = sqrt(v0.x()*v0.x() + v0.y()*v0.y() + v0.z()*v0.z());
  const double r1 = sqrt(v1.x()*v1.x() + v1.y()*v1.y() + v1.z()*v1.z());
  const double dr = r1-r0;
  double s = lcmin;
  int n = (int)ceil(sqrt(dx*dx+dy*dy+dz*dz)/(0.5*s));
  for (int i = 0; i < n-1; ++i) {
    double t = (i+1)/(double)n;
    double xyz[3] = {v0.x()+t*dx,v0.y()+t*dy,v0.z()+t*dz};
    radiusScale(r0+t*dr, xyz);
    vertices.push_back(new Vertex(0,xyz[0],xyz[1],xyz[2],-1,tag));
  }
  double t = (0.1)/(double)n;
  if(withShiftedStart){
    double xyz[3] = {v0.x()+t*dx,v0.y()+t*dy,v0.z()+t*dz};
    radiusScale(r0+t*dr, xyz);
    vertices.push_back(new Vertex(0,xyz[0],xyz[1],xyz[2],-1,tag));
  }
  t = 1-t;
  if(withShiftedEnd){
    double xyz[3] = {v0.x()+t*dx,v0.y()+t*dy,v0.z()+t*dz};
    radiusScale(r0+t*dr, xyz);
    vertices.push_back(new Vertex(0,xyz[0],xyz[1],xyz[2],-1,tag));
  }
}


//see http://www.oocities.org/geoff_wass/dBASE/GaryWhite/dBASE/FAQ/qformt.htm#A
std::vector<int> readPhysicalsFromDBF(const char *shpfilename) {
  std::vector<int> physicals;
  std::string shpfname(shpfilename);
  if (shpfname.size() < 4) return physicals;
  std::string dbffname = shpfname.substr(0,shpfname.size()-4) + ".dbf";
  std::ifstream dbffile(dbffname.c_str(), std::ios::binary|std::ios::in);
  if (!dbffile.is_open())
    return physicals;
  char dbh[32];
  try {
    if (!dbffile.read(dbh,32)) throw 0;
    unsigned char format = dbh[0];
    if (format != 3) throw 0;
    int nrecord = *(int*)(dbh+4);
    physicals.resize(nrecord);
    int recordsize = *(int*)(dbh+10);
    int i = 0;
    int cur = 1;
    int pStart=-1, pLength;
    while(true) {
      if (!dbffile.read(dbh,32)) throw 1;
      std::string fieldname(dbh + 0);
      char type = dbh[11];
      unsigned char length = dbh[16];
      if ((fieldname == "entity" || fieldname == "physical") && type == 'N'){
        pStart = cur;
        pLength = length;
      }
      cur += length;
      char isEnd;
      if (!dbffile.read((char*)&isEnd,1)) throw 2;
      if (isEnd == 0xd || dbffile.eof())
        break;
      dbffile.unget();
    }
    if (pStart == -1) {
      printf("dbf file has no \"entity\" field\n");
      throw 3;
    }
    std::vector<char> val(recordsize);
    for (int i = 0; i < nrecord; ++i) {
      char present;
      if (!dbffile.read(&val[0],recordsize)) throw 4;
      if (val[0] == '*') continue;
      if (val[0] != 0x20) ERROR;
      try {
        physicals[i] = std::stoi(std::string(&val[pStart], &val[pStart+pLength]));
      }
      catch(int e) {
        physicals[i] = 0;
      }
    }
  }
  catch (int e) {
    printf("invalid dbf file (error %i)\n", e);
  }
  return physicals;
}

std::vector<int> readSHX(const char *shpfilename) {
  std::vector<int> starts;
  std::string shpfname(shpfilename);
  if (shpfname.size() < 4) return starts;
  std::string shxfname = shpfname.substr(0,shpfname.size()-4) + ".shx";
  std::ifstream file(shxfname.c_str(), std::ios::binary|std::ios::in);
  if (!file.is_open())
    return starts;
  int fileCode, fileLength, fileVersion, shapeType;
  if(!file.read((char*)&fileCode, 4)) ERROR;
  swapIntByte(fileCode);
  if(fileCode != 9994) ERROR;
  file.seekg(20, std::ios_base::cur);
  if(!file.read((char*)&fileLength, 4)) ERROR;
  swapIntByte(fileLength);
  fileLength*=2;
  if(!file.read((char*)&fileVersion, 4)) ERROR;
  if(!file.read((char*)&shapeType, 4)) ERROR;
  if(shapeType != 3 && shapeType !=5 && shapeType !=13) ERROR; // PolyLine or Polygone
  int lread = 100;
  file.seekg(lread, std::ios_base::beg);
  while(lread < fileLength) {
    int start, length;
    if(!file.read((char*)&start, 4)) ERROR;
    swapIntByte(start);
    starts.push_back(start);
    if(!file.read((char*)&length, 4)) ERROR;
    //swapIntByte(length);
    lread += 8;
  }
  return starts;
}

static std::vector<Vertex*> readSHP(const char *filename, double lcmin, std::string ellps) {
  projPJ pjglobal = pj_init_plus(("+proj=geocent "+ellps).c_str());
  projPJ pjll = pj_init_plus(("+proj=lonlat "+ellps).c_str());
  if(!pjglobal)
  {
    printf("Invalid ellipsoid : %s\n", ellps.c_str());
    exit(1);
  }
  std::vector<Vertex*> vertices;
  std::vector<int> physicals = readPhysicalsFromDBF(filename);
  std::vector<int> starts = readSHX(filename);
  std::ifstream file(filename,std::ios::binary|std::ios::in);
  int fileCode, fileLength, fileVersion, shapeType;
  if(!file.read((char*)&fileCode, 4)) ERROR;
  swapIntByte(fileCode);
  if(fileCode != 9994) ERROR;
  file.seekg(20, std::ios_base::cur);
  if(!file.read((char*)&fileLength, 4)) ERROR;
  swapIntByte(fileLength);
  fileLength*=2;
  if(!file.read((char*)&fileVersion, 4)) ERROR;
  if(!file.read((char*)&shapeType, 4)) ERROR;
  if(shapeType != 3 && shapeType !=5 && shapeType !=13) ERROR; // PolyLine or Polygone
  int lread = 100;
  int irecord = 0;
  while(lread < fileLength) {
    if(!starts.empty()){
      if(irecord > starts.size())
        break;
      lread = starts[irecord]*2;
    }
    file.seekg(lread, std::ios_base::beg);
    int id, length, type;
    double bbox[4];
    if(!file.read((char*)&id, 4)) ERROR;
    swapIntByte(id);
    int physical = physicals.size() > id-1 ? physicals[id-1] : 0;
    if(!file.read((char*)&length, 4)) ERROR;
    swapIntByte(length);
    if(!file.read((char*)&type, 4)) ERROR;
    if(type != 3 && type !=5 && type !=13) ERROR;
    file.seekg(32,std::ios_base::cur);
    int npart, npoint;
    if(!file.read((char*)&npart, 4)) ERROR;
    if(!file.read((char*)&npoint, 4)) ERROR;
    lread += 2*length + 8;
    std::vector<int> part(npart+1);
    if(!file.read((char*)&part[0], 4*npart)) ERROR;
    part.back() = npoint;
    for(size_t i = 0; i < part.size()-1; ++i) {
      std::vector<double> points(2*(part[i+1]-part[i]));
      if(!file.read((char*)&points[0], 8*points.size())) ERROR;
      size_t cur = vertices.size();
      std::vector<Vertex*> lv;
      for (size_t j = 0; j < (points.size()/2); ++j) {
        double xyz[3] = {points[j*2]*M_PI/180, points[j*2+1]*M_PI/180, 0};
        pj_transform(pjll, pjglobal, 1, 1, &xyz[0], &xyz[1], &xyz[2]);
        lv.push_back(new Vertex(0, xyz[0], xyz[1], xyz[2], -1., physical));
      }
      for (size_t j = 0; j < lv.size()-1; ++j) {
        saturateEdge(vertices, *lv[j], *lv[j+1], lcmin, physical, j==0, j==lv.size()-2);
        if (j>0){
          vertices.push_back(lv[j]);
        }
      }
    }
    irecord++;
  }
  return vertices;
}

extern "C" {
  void read_shp(const char *filename, double lcmin, const char *ellps,int *n, double **x, int **tag);
  void libcfree(void*p);
}

void read_shp(const char *filename, double lcmin, const char *ellps,int *n, double **x, int **tag) {
  auto v = readSHP(filename,lcmin,ellps);
  *n = v.size();
  *x = (double*)malloc(sizeof(double)*v.size()*3);
  *tag = (int*)malloc(sizeof(int)*v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    (*x)[i*3+0] = v[i]->x();
    (*x)[i*3+1] = v[i]->y();
    (*x)[i*3+2] = v[i]->z();
    (*tag)[i] = v[i]->tag();
  }
}

void libcfree(void *p) {
  free(p);
}
