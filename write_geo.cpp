
void writeGeo(faceContainer &mesh, const char *filename, double lc, std::map<int,std::string> boundary_names, projectorCB proj, void *projdata, int NT)
{
  std::ofstream f(filename);
  std::map<std::string,std::vector<int> > physicals;
  f << "Mesh.LcIntegrationPrecision = 1e-4;\n";
  f << "LC="<<lc<<";\n";
  size_t ip = 1;
  if (!proj){
    f <<"Point("<<ip++<<")={0,0,0};\n";
    f <<"Point("<<ip++<<")={0,0,-1};\n";
    f << "PolarSphere(1)={1,2};\n";
  }
  size_t il = 1;
  for (int thread=0;thread<NT;thread++){
    for (int i=0;i<mesh.size(thread);i++){
      mesh(thread,i)->_modified = false;
    }
  }
  size_t ill = 1;
  for (int thread=0;thread<NT;thread++){ for (int i=0;i<mesh.size(thread);i++){
      Face &T = *mesh(thread,i);
      if(T._color != 0 || !T.V[0] || T._modified) continue;
      for (int j=0; j<3; ++j) {
        if (T.F[j]->_color == 0|| T._modified) continue;
        int firstp = ip;
        int firstpll = ip;
        int firstl = il;
        Vertex *firstV = T.V[j], *curV = T.V[j];
        Face *curF = &T;
        int iV = getid3(curF->V, curV);
        int oldtag = -1;
        do {
          curF->_modified = true;
          if (proj) {
            double x[3] = {curV->x(), curV->y(), curV->z()};
            proj(x,projdata);
            f << "Point(" << ip++ <<")={"<<x[0] <<","<<x[1]<<","<<x[2]<<",LC};\n";
          }
          else {
            f << "Point(" << ip++ <<")={"<<-curV->x()/(1+curV->z()) <<","<<-curV->y()/(1+curV->z())<<","<<0<<"};\n";
          }
          int tag = ((Vertex*)curV)->tag();
          if (oldtag != -1 && oldtag != tag) {
            physicals[boundary_names[oldtag]].push_back(il);
            f << "Spline(" << il++ << ")={"<< firstp <<":"<<ip-1<<"};\n";
            firstp = ip-1;
          }
          oldtag = tag;

          iV = (iV+1)%3;
          curV = curF->V[iV];
          while(curF->F[iV]->_color == 0) {
            curF = curF->F[iV];
            iV = getid3(curF->V,curV);
          }
        }while(curV != firstV);
        physicals[boundary_names[oldtag]].push_back(il);
        f << "Spline(" << il++ << ")={"<< firstp <<":"<<ip-1<<","<<firstpll<<"};\n";
        f << "Line Loop(" << ill++ << ")={"<< firstl << ":" << il-1 << "};\n";
        il++;
      }
    }
  }
  f << "Physical Surface(\"domain\") = {1};\n";
  for(auto p: physicals) {
    f << "Physical Line(\"" << p.first << "\") = {";
    for(size_t i = 0; i < p.second.size()-1; ++i) {
      f << p.second[i] << ",";
    }
    f << p.second.back() << "};\n";

  }
  f << "Plane Surface(1) = {1:"<<ill-1<<"};\n";
  f << "Mesh.CharacteristicLengthExtendFromBoundary = 0;\n";
  f << "Mesh.CharacteristicLengthFromPoints = 0;\n";
}
