#ifndef VERTEX_H_
#define VERTEX_H_
#include <climits>
const int MAX_NUM_THREADS_ = 16;
const size_t OFFSET_NUM_ = ULONG_MAX / MAX_NUM_THREADS_;
struct Face;
struct Vertex {
private :
  double _x[3];
  double _lc;
  unsigned char _thread;
  Face *_f;
  int _entity;
public :
  char _color;
  static double sizeField(double x, double y, double z);
  inline double x() const {return _x[0];}
  inline double y() const {return _x[1];}
  inline double z() const {return _x[2];}
  inline double lc() const {return _lc;}
  inline double &x() {return _x[0];}
  inline double &y() {return _x[1];}
  inline double &z() {return _x[2];}
  inline double &lc() {return _lc;}
  inline void  setF (Face*f) {_f=f;}
  inline Face* getF () const {return _f;}
  inline operator double *() { return _x; }
  int tag() {return _entity;}
  void setTag(int tag) {_entity = tag;}
  Vertex (int thread, double X, double Y, double Z, double lc, int tag=0) : _thread(thread), _f(0){
    _x[0]=X; 
    _x[1]=Y;
    _x[2]=Z;
    _lc=lc;
    _color = 0;
    _entity = tag;
  }
  //inline SPoint3 point() const { return SPoint3(x(), y(), z()); }
};
#endif
