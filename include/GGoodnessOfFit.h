#ifndef GGOODNESSOFFIT
#define GGOODNESSOFFIT

#include <vector>
#include <string>
namespace Gamapola
{
    class GGoodnessOfFit
    {
    public:
        virtual ~GGoodnessOfFit();
        virtual void GMakeTest() = 0;
        virtual void GSetTruthEvent(const std::vector<double>& phaseSpace) = 0;
        virtual void GSetGeneratedData(const std::vector<double>& modelPars, const int& nEvts, const int& charge, const std::string& cut) = 0;
        virtual const std::vector<double>& GGetpValues() const = 0;
//         virtual void GSetModelParsErrors(const std::vector<double>& modelParsErrors) = 0;
    };
}

#endif
