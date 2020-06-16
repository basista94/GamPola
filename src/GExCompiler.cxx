#include "GExCompiler.h"

namespace Gamapola{
  GExCompiler::GExCompiler():
  fArg(0)
  {
    std::cout << "GExCompiler constructor calling. . ." << std::endl;
  }
  GExCompiler::~GExCompiler()
  {
//     std::cout << "GExCompiler destructor calling. . ." << std::endl;
  }
  void GExCompiler::GCreateSrcFile(std::string& filename, std::ofstream& ofs)
  {
    if (filename.empty()) 
    {
      const char* filename_pattern = "./GiNaCXXXXXX";
      char* new_filename = new char[strlen(filename_pattern)+1];
      strcpy(new_filename, filename_pattern);
      if (!mkstemp(new_filename)) 
      {
        delete[] new_filename;
        throw std::runtime_error("mktemp failed");
      }
      filename = std::string(new_filename);
      ofs.open(new_filename, std::ios::out);
      delete[] new_filename;
     } 
     else 
      ofs.open(filename.c_str(), std::ios::out);
                
     if (!ofs) 
      throw std::runtime_error("could not create source code file for compilation");

     ofs << "#include <stddef.h> " << std::endl;
     ofs << "#include <stdlib.h> " << std::endl;
     ofs << "#include <math.h> " << std::endl;
     ofs << "#include <iostream> " << std::endl;
     ofs << "#include <complex> " << std::endl;
     ofs << "#ifdef __cplusplus " << std::endl;
     ofs << "extern \"C\" { " << std::endl;
     ofs << "#endif " << std::endl;
     ofs << std::endl;
  }
  void GExCompiler::GCompileSrcFile(const std::string filename, bool clean_up)
  {
    std::string strcompile = "../additionalFiles/ginac-excompiler2 " + filename;
    if (system(strcompile.c_str())) 
      throw std::runtime_error("excompiler::compile_src_file: error compiling source file!");
    if (clean_up) 
      remove(filename.c_str());
  }
  void* GExCompiler::GLinkSoFile(const std::string filename, bool clean_up)
  {
    void* module = NULL;
    module = dlopen(filename.c_str(), RTLD_NOW);
    if (module == NULL)
      throw std::runtime_error("excompiler::link_so_file: could not open compiled module!");
//     add_opened_module(module, filename, clean_up);
    if(dlsym(module, "compiled_ex"))
      std::cout << "Opening module " << filename << std::endl;
    return dlsym(module, "compiled_ex");
  }
  void GExCompiler::GCompileEx(const ex& expr, const symbol& sym, GFUNCP_1P& fp, const std::string filename)
  {
        symbol x("x");
        ex expr_with_x = expr.subs(lst(sym==x));

        std::ofstream ofs;
        std::string unique_filename = filename;
        GCreateSrcFile(unique_filename, ofs);

        ofs << "double compiled_ex(double x)" << std::endl;
        ofs << "{" << std::endl;
        ofs << "double res = ";
        expr_with_x.print(GiNaC::print_csrc_double(ofs));
        ofs << ";" << std::endl;
        ofs << "return(res); " << std::endl;
        ofs << "}" << std::endl;

        ofs << "#ifdef __cplusplus " << std::endl;
        ofs << " } " << std::endl;
        ofs << "#endif " << std::endl;
        
        ofs.close();

        GCompileSrcFile(unique_filename, filename.empty());
        // This is not standard compliant! ... no conversion between
        // pointer-to-functions and pointer-to-objects ...
        fp = (GFUNCP_1P) GLinkSoFile("./"+unique_filename+".so", filename.empty());
  }
  void GExCompiler::GCompileEx(const lst& exprs, const lst& syms, GFUNCP_CUBA& fp, const std::string filename)
  {
        lst replacements;
        for (std::size_t count=0; count<syms.nops(); ++count) {
                std::ostringstream s;
                s << "a[" << count << "]";
                replacements.append(syms.op(count) == symbol(s.str()));
        }

        std::vector<ex> expr_with_cname;
        for (std::size_t count=0; count<exprs.nops(); ++count) {
                expr_with_cname.push_back(exprs.op(count).subs(replacements));
        }

        std::ofstream ofs;
        std::string unique_filename = filename;
        GCreateSrcFile(unique_filename, ofs);

        ofs << "void compiled_ex(std::complex<double>* a, std::complex<double>* f)" << std::endl;
        ofs << "{" << std::endl;
//         std::cout << exprs.nops() << std::endl;
//         ofs << "        std::complex<double> sinTh = sin(a[2]);" << std::endl;
//         ofs << "        std::complex<double> cosTh = cos(a[2]);" << std::endl;
        for (std::size_t count=0; count<exprs.nops(); ++count) {
//           std::cout << exprs.op(count) << std::endl;
                ofs << "        f[" << count << "] = ";
                expr_with_cname[count].print(GiNaC::print_csrc_double(ofs));
                ofs << ";" << std::endl;
//                 ofs << "        std::cout << f[" << count << "] << \"   \" <<  " << count << " << std::endl;" << std::endl;
        }
        ofs << "}" << std::endl;
        
        ofs << "#ifdef __cplusplus " << std::endl;
        ofs << " } " << std::endl;
        ofs << "#endif " << std::endl;
        ofs.close();
        GCompileSrcFile(unique_filename, filename.empty());
//         std::cout << unique_filename << '\n';
        fp = (GFUNCP_CUBA) GLinkSoFile("./"+unique_filename+".so", filename.empty());
//         if(fp)
//           std::cout << "AAAAAAAAAAAAAAAAAAAAa" << std::endl;
  }
}