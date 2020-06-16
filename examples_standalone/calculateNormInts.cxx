#include "GFitter.h"

using namespace Gamapola;

int main()
{
    auto f = std::make_shared<GFitter>(0, "dummy.csv");
    f->GGenerateNormalizationIntegrals("norm_neutral_5k_cut1_1-1.6.txt", "1 < sKpipi < 1.6");
    return 0;
}

