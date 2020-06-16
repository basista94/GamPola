#ifndef GEXCOMPILER_H
#define GEXCOMPILER_H

#include <iostream>
#include "ginac/ginac.h"
#include <dlfcn.h>
#include <fstream>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>

using namespace GiNaC;
namespace Gamapola
{
  typedef double (*GFUNCP_1P) (double);
  typedef double (*GFUNCP_2P) (double, double);
  typedef void (*GFUNCP_CUBA) (std::complex<double>* , std::complex<double>* );
  
  class GExCompiler{
  public:
    GExCompiler();
    virtual ~GExCompiler();
    void GCompileEx(const ex& expr, const symbol& sym, GFUNCP_1P& fp, const std::string filename = "");
    void GCompileEx(const ex& expr, const symbol& sym1, const symbol& sym2, GFUNCP_2P& fp, const std::string filename = "");
    void GCompileEx(const lst& exprs, const lst& syms, GFUNCP_CUBA& fp, const std::string filename = "");
    void GLinkEx(const std::string filename, GFUNCP_1P& fp);
    void GLinkEx(const std::string filename, GFUNCP_2P& fp);
    void GLinkEx(const std::string filename, GFUNCP_CUBA& fp);
    void GUnlinkEx(const std::string filename);
    void GCreateSrcFile(std::string& filename, std::ofstream& ofs);
    void GCompileSrcFile(const std::string filename, bool clean_up);
    void* GLinkSoFile(const std::string filename, bool clean_up);
  private:
    int fArg;
  };
}
#endif