#include "GBF2Ana.h"
#include <memory>
#include <unordered_map>

using namespace Gamapola;

int main()
{
    std::unordered_map<std::string, std::tuple<int, double, double>> dict; 
    dict["Mkpipi"] = std::tuple<int, double, double>{100, 1.1, 1.9};
    dict["Mpipi"] = std::tuple<int, double, double>{44, 0.225, 1.525};
    dict["Mkpi2"] = std::tuple<int, double, double>{37, 0.58, 1.66};
    dict["cosTheta"] = std::tuple<int, double, double>{20, -1, 1};
    
    auto&& gbf2 = std::make_shared<GBF2Ana>(1, "../additionalFiles/model4BF2.fit", "norm_charged_5k_cut1_1-1_9.txt");
    gbf2->GHistogramize(dict, "../additionalFiles/model.ge", "histogram.csv", 100000);
    gbf2->GSetData("../lhcb_mkpipi.csv", "Mkpipi");
    gbf2->GSetData("../lhcb_mpipi.csv", "Mpipi");
    gbf2->GSetData("../lhcb_mkpi.csv", "Mkpi2");
    gbf2->GFit("../additionalFiles/model4BF2.fit");
    return 0;
}
