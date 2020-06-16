#include "EvtGenModels/EvtGSimpleRandomEngine.hh"

double GSimpleRandomEngine::random(){
  
  _next=_next*1103515245+123345;
  unsigned temp=(unsigned)(_next/65536) % 32768;
  
  return ( temp + 1.0 ) / 32769.0;
}



