#include "GFitter.h"

using namespace Gamapola;

int main()
{  
    auto g1 = std::make_shared<GGenerator>(0);
    g1->GSetCouplings("../additionalFiles/model.ge");
    g1->GGenerate(1000);
    g1->GWriteToFile("neutral_1000.root");
    
    auto f = std::make_shared<GFitter>(0, "wrong_model_stat_fit_1.6.csv");
    f->GSetModelParameters("../additionalFiles/model.fit");    
    f->GSetNormalizationIntegrals("norm_neutral_5k_cut1_1-1.6.txt");
    f->GSetData("neutral_1000.root", "1 < sKpipi < 1.6");
    f->GFit(false);
    
    return 0;
}


