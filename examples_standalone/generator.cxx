#include "GGenerator.h"

using namespace Gamapola;

int main()
{
    auto g1 = std::make_shared<GGenerator>(0);
    g1->GSetCouplings("../additionalFiles/model.ge");
    g1->GGenerate(10000);
    g1->GWriteToFile("neutral_100k.root");
    return 0;
}
