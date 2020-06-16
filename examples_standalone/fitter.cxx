#include "GFitter.h"

using namespace Gamapola;

int main()
{  
        auto f = std::make_shared<GFitter>(0, "dummy.csv");
        f->GSetModelParameters("../additionalFiles/model.fit");    
        f->GSetNormalizationIntegrals("norm_neutral_5k_cut1_1-1.6.txt");
        f->GSetData("neutral_100k.root", "1 < sKpipi < 1.6");
        for(auto i = 0; i < 10; ++i)
            f->GFit(true);
        return 0;
}


