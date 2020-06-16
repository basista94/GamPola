#ifndef GWRITETOFILE_H
#define GWRITETOFILE_H

#include <iostream>

namespace Gamapola
{
  class GWriteToFile{
  public:
    GWriteToFile();
    virtual ~GWriteToFile();
  private:
    int fArg;
  };
}
#endif