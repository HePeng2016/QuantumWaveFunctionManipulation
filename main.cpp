 gsl_Type Context;
 char * buffer;
 buffer = (char *)malloc(sizeof(char)*1024);
 InitInvariantContext(&Context);
 std::vector< Coordinate > AtomN_o;
 std::vector<double> MO_o;
 InvariantDirection(bs,MOs[5],AtomN,MO_o,AtomN_o,&Context);
 FreeInvariantContext(&Context);
